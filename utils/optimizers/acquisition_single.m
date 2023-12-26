function [solution_promising,improvement,paras] = acquisition_single(database,lb,ub,name_saea,paras)

switch name_saea

    %**************Bayesian Optimization with Lower Confidence Bound**************%
    case 'BO-LCB'
        % build the surrogate model using GPR
        flag = 0;
        while flag == 0
            try
                func_surrogate = surrogate_model(database,'gpr');
                flag = 1;
            catch
                continue;
            end
        end
        % perform the evolutionary search on the surrogate model
        configs.database = database;
        configs.popsize = 50;
        configs.query_max = 50;
        configs.ini_strategy = 'elite+random';
        configs.acquisition = 'lcb';
        configs.search_engine = 'ea';
        configs.selector = 'roulette_wheel';
        [population,objs_acq] = acquisition_template(func_surrogate,lb,ub,configs);
        % estimate the internal improvement
        [~,idx] = min(objs_acq);
        solution_promising = population(idx,:);
        improvement = improvement_internal(func_surrogate,configs.acquisition,database(:,1:end-1),solution_promising,size(database,1));
        %**************Bayesian Optimization with Lower Confidence Bound**************%



        %********************Three-Level Radial Basis Function Method********************%
        % Reference: Li, G., Zhang, Q., Lin, Q., & Gao, W. (2021). A three-level radial basis function
        % method for expensive optimization. IEEE Transactions on Cybernetics, 52(7), 5720-5731.
    case 'TLRBF'
        % initialization
        dim = size(database,2)-1;
        no_initials = 50;
        state = mod(size(database,1)-no_initials,3); % 0: global; 1: medium; 2: local.
        model = 'rbf';
        % optimization
        switch(state)
            case 0 % global search
                alpha = 0.4;
                m = 200*dim;
                func_surrogate_global = surrogate_model(database,model);
                solutions_global = repmat(lb,m,1)+(repmat(ub,m,1)-repmat(lb,m,1)).*rand(m,dim);
                dist=pdist2(solutions_global,database(:,1:end-1));
                mindist=min(dist,[],2);
                deltag=alpha*max(mindist);
                index=find(mindist<=deltag);
                if ~isempty(index)
                    solutions_global(index,:)=[];
                end
                objs_pre_global = func_surrogate_global(solutions_global);
                [~,idx] = min(objs_pre_global);
                solution_promising = solutions_global(idx,:);
                improvement = min(database(:,end))-min(objs_pre_global);
            case 1 % subregion search
                N = size(database,1);
                L1 = 100;
                L2 = 80;
                if N <= L1
                    X_subregion = database(:,1:end-1);
                    Y_subregion = database(:,end);
                    func_surrogate_subregion = surrogate_model(database,model);
                else
                    X = database(:,1:end-1);
                    Y = database(:,end);
                    no_clusters = 1+ceil((N-L1)/L2);
                    X_normalized = mapminmax(X',0,1)';
                    [~,U,~] = fcm_tlrbf(X_normalized,no_clusters);
                    U = U';
                    [~,idx] = sort(U,'descend');
                    X_cluster=[];
                    mper=[];
                    Y_cluster = [];
                    for k=1:no_clusters
                        X_cluster{k}=X(idx(1:L1,k),:);
                        Y_cluster{k}=Y(idx(1:L1,k),1);
                        mper(k)=mean(Y_cluster{k});
                        func_surrogate{k} = surrogate_model([X_cluster{k},Y_cluster{k}],model);
                    end
                    [~,idx]=sort(mper,'descend');
                    pro=idx/no_clusters;
                    sid=randi([1,no_clusters],1);
                    while rand>pro(sid)
                        sid=randi([1,no_clusters],1);
                    end
                    X_subregion = X_cluster{sid};
                    Y_subregion = Y_cluster{sid};
                    func_surrogate_subregion = func_surrogate{sid};
                end
                % perform the evolutionary search on the subregion surrogate model
                configs.database = [X_subregion Y_subregion];
                configs.popsize = 50;
                configs.query_max = 50;
                configs.ini_strategy = 'elite+random';
                configs.acquisition = 'plain';
                configs.search_engine = 'ea';
                configs.selector = 'roulette_wheel';
                lb_subregion = min(X_subregion);
                ub_subregion = max(X_subregion);
                [population,objs_acq] = acquisition_template(func_surrogate_subregion,lb_subregion,ub_subregion,configs);
                [~,idx] = min(objs_acq);
                solution_promising = population(idx,:);
                improvement = min(Y_subregion)-func_surrogate_subregion(solution_promising);
            case 2 % local search
                k = min(2*dim,size(database,1)-1);
                X = database(:,1:end-1);
                Y = database(:,end);
                [~,idx_min] = min(Y);
                dist = pdist2(X,X(idx_min,:));
                [~,idx] = sort(dist,'ascend');
                sid = idx(1:k);
                X_local = X(sid,:);
                Y_local = Y(sid,1);
                lb_local = min(X_local);
                ub_local = max(X_local);
                func_surrogate_local = surrogate_model([X_local,Y_local],model);
                % perform the evolutionary search on the subregion surrogate model
                configs.database = [X_local Y_local];
                configs.popsize = 50;
                configs.query_max = 50;
                configs.ini_strategy = 'elite+random';
                configs.acquisition = 'plain';
                configs.search_engine = 'ea';
                configs.selector = 'roulette_wheel';
                [population,objs_acq] = acquisition_template(func_surrogate_local,lb_local,ub_local,configs);
                [~,idx] = min(objs_acq);
                solution_promising = population(idx,:);
                improvement = min(Y_local)-func_surrogate_local(solution_promising);
        end
        %********************Three-Level Radial Basis Function Method********************%



        %*************Global and Local Surrogate-Assisted Differential Evolution*************%
        % Reference: Wang, W., Liu, H. L., & Tan, K. C. (2022). A surrogate-assisted differential
        % evolution algorithm for high-dimensional expensive optimization problems. IEEE
        % Transactions on Cybernetics, 53(4), 2685-2697.
    case 'GL-SADE'
        % initialization
        dim = size(database,2)-1;
        if database(end,end) < min(database(1:end-1,end))
            state = 1;
        else
            state = 0;
        end
        model_global = 'rbf';
        model_local = 'gpr';
        switch(state)
            case 0 % global search
                func_surrogate_global = surrogate_model(database,model_global);
                % perform the evolutionary search on the global surrogate model
                configs.database = database;
                configs.popsize = 50;
                configs.query_max = 50;
                configs.ini_strategy = 'elite';
                configs.acquisition = 'plain';
                configs.search_engine = 'de-b';
                configs.selector = 'comparison';
                [population,objs_acq] = acquisition_template(func_surrogate_global,lb,ub,configs);
                % estimate the internal improvement
                [~,idx] = min(objs_acq);
                solution_promising = population(idx,:);
                improvement = improvement_internal(func_surrogate_global,configs.acquisition,database(:,1:end-1),solution_promising,size(database,1));
            case 1 % local search
                no_points_local = min([size(database,1),2*dim]);
                X_local = database(end-no_points_local+1:end,1:end-1);
                Y_local = database(end-no_points_local+1:end,end);
                func_surrogate_local = surrogate_model([X_local,Y_local],model_local);
                % perform the evolutionary prescreening on the local surrogate model
                configs.database = database;
                configs.popsize = 50;
                configs.query_max = 1; % prescreening
                configs.ini_strategy = 'elite';
                configs.acquisition = 'lcb-decay';
                configs.search_engine = 'de-b';
                configs.selector = 'comparison';
                [population,objs_acq] = acquisition_template(func_surrogate_local,lb,ub,configs);
                % estimate the internal improvement
                [~,idx] = min(objs_acq);
                solution_promising = population(idx,:);
                improvement = improvement_internal(func_surrogate_local,configs.acquisition,X_local,solution_promising,size(database,1));
        end
        %*************Global and Local Surrogate-Assisted Differential Evolution*************%



        %************************Multi-Evolutionary Sampling Strategy************************%
        % Reference: Yu, F., Gong, W., & Zhen, H. (2022). A data-driven evolutionary algorithm
        % with multi-evolutionary sampling strategy for expensive optimization.
        % Knowledge-Based Systems, 242, 108436.
    case 'DDEA-MESS'
        % initialization
        strategy_id = mess(size(database,1),500);
        if strategy_id == 1 % global search
            m = min([size(database,1),300]);
            func_surrogate_global = surrogate_model(database(1:m,:),'rbf');
            % perform the evolutionary prescreening on the global surrogate model
            configs.database = database;
            configs.popsize = 50;
            configs.query_max = 1; % prescreening
            configs.ini_strategy = 'elite';
            configs.acquisition = 'plain';
            configs.search_engine = 'de-r';
            configs.selector = 'comparison';
            [population,objs_acq] = acquisition_template(func_surrogate_global,lb,ub,configs);
            % estimate the internal improvement
            [~,idx] = min(objs_acq);
            solution_promising = population(idx,:);
            improvement = min(database(:,end))-func_surrogate_global(solution_promising);
        elseif strategy_id == 2 % local search
            tao = min([size(database,2)-1+25,size(database,1)]);
            [~,idx] = sort(database(:,end));
            X_local = database(idx(1:tao),1:end-1);
            Y_local = database(idx(1:tao),end);
            lb_local = min(X_local);
            ub_local = max(X_local);
            func_surrogate_local = surrogate_model([X_local, Y_local],'rbf');
            % perform the evolutionary search on the local surrogate model
            configs.database = database;
            configs.popsize = 50;
            configs.query_max = 10;
            configs.ini_strategy = 'elite';
            configs.acquisition = 'plain';
            configs.search_engine = 'de-b';
            configs.selector = 'comparison';
            [population,objs_acq] = acquisition_template(func_surrogate_local,lb_local,ub_local,configs);
            % estimate the internal improvement
            [~,idx] = min(objs_acq);
            solution_promising = population(idx,:);
            improvement = min(Y_local)-func_surrogate_local(solution_promising);
        else % interior point search
            m = min([size(database,1),5*(size(database,2)-1)]);
            X = database(:,1:end-1);
            Y = database(:,end);
            [~,idx_min] = min(Y);
            dist = pdist2(X,X(idx_min,:));
            [~,idx] = sort(dist,'ascend');
            sid = idx(1:m);
            X_trs = X(sid,:);
            Y_trs = Y(sid,1);
            func_surrogate_tr = surrogate_model([X_trs,Y_trs],'rbf');
            options = optimoptions('fmincon','Algorithm','interior-point','MaxIterations',20,'Display','off');
            solution_promising= fmincon(func_surrogate_tr,X(idx_min,:),[],[],[],[],min(X_trs),max(X_trs),[],options);
            improvement = min(Y_trs)-func_surrogate_tr(solution_promising);
        end
        %************************Multi-Evolutionary Sampling Strategy************************%



        %******************Lipschitz Surrogate-assisted Differential Evolution******************%
        % Reference: Kůdela, J., & Matoušek, R. (2023). Combining Lipschitz and RBF surrogate
        % models for high-dimensional computationally expensive problems. Information
        % Sciences, 619, 457-477.
    case 'LSADE'
        % initialization
        no_initials = 50;
        state = mod(size(database,1)-no_initials,3); % static LSADE: 0: RBF; 1: Lipschitz; 2: Local.
        switch(state)
            case 0 % RBF condition
                func_surrogate_rbf = surrogate_model(database,'rbf');
                popsize = 50;
                acquisition = 'plain';
                ini_strategy = 'best+randperm';
                search_engine = 'de-b';
                [population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,[],[]); % initialize the population
                population_offspring = pop_offspring(population,objs_acq,lb,ub,search_engine); % generate the offspring population
                objs_acq_offspring = pop_acquisition(func_surrogate_rbf,population_offspring,acquisition,size(database,1)); % acquisite the objective values of the offspring
                % estimate the internal improvement
                [~,idx] = min(objs_acq_offspring);
                solution_promising = population_offspring(idx,:);
                improvement = improvement_internal(func_surrogate_rbf,acquisition,database(:,1:end-1),solution_promising,size(database,1));
            case 1 % Lipschitz condition
                k = k_est(database(:,end),database(:,1:end-1)');
                func_surrogate_lip = @(x)f_lip(database(:,end),database(:,1:end-1)',k,x');
                popsize = 50;
                acquisition = 'plain';
                ini_strategy = 'best+randperm';
                search_engine = 'de-b';
                [population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,[],[]); % initialize the population
                population_offspring = pop_offspring(population,objs_acq,lb,ub,search_engine); % generate the offspring population
                objs_acq_offspring = pop_acquisition(func_surrogate_lip,population_offspring,acquisition,size(database,1)); % acquisite the objective values of the offspring
                % estimate the internal improvement
                [~,idx] = min(objs_acq_offspring);
                solution_promising = population_offspring(idx,:);
                improvement = improvement_internal(func_surrogate_lip,acquisition,database(:,1:end-1),solution_promising,size(database,1));
            case 2 % Local optimization condition
                no_points_local = min([size(database,1),3*(size(database,2)-1)]);
                [~,idx] = sort(database(:,end));
                X_local = database(idx(1:no_points_local),1:end-1);
                Y_local = database(idx(1:no_points_local),end);
                lb_local = min(X_local);
                ub_local = max(X_local);
                func_surrogate = surrogate_model([X_local,Y_local],'rbf');
                options = optimoptions("fmincon",'Algorithm','sqp','maxiter',20,'Display','off');
                x_initial = lb_local + (ub_local-lb_local).*rand(1,size(database,2)-1);
                [solution_promising,obj_local] = fmincon(func_surrogate,x_initial,[],[],[],[],lb_local,ub_local,[],options);
                % estimate the internal improvement
                improvement = min(Y_local)-obj_local;
        end
        %******************Lipschitz Surrogate-assisted Differential Evolution******************%



        %********************Model and Infill Criterion Auto-Configuration********************%
        % Reference: Xie, L., Li, G., Wang, Z., Cui, L., & Gong, M. (2023). Surrogate-assisted
        % evolutionary algorithm with model and infill criterion auto-configuration. IEEE
        % Transactions on Evolutionary Computation.
    case 'AutoSAEA'
        no_initials = 50;
        no_arms = 8; % Number of combinatorial arms
        id = size(database,1)-no_initials+1;
        popsize = 50;
        ini_strategy = 'elite';
        search_engine = 'de-b';
        [population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,[],[]); % initialize the population
        population_offspring = pop_offspring(population,objs_acq,lb,ub,search_engine); % generate the offspring population
        if id <= no_arms
            if id == 1 % {RBF, prescreening}
                [solution_promising, improvement,reward_rp] = rbf_pre_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                paras.save_rp = [paras.save_rp, reward_rp];
                paras.rbf_model = [paras.rbf_model;reward_rp];
            elseif id == 2 % {GP, LCB}
                [solution_promising, improvement,reward_gl] = gp_lcb_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                paras.save_gl = [paras.save_gl, reward_gl];
                paras.gp_model = [paras.gp_model;reward_gl];
            elseif id == 3 % {RBF, local search}
                [solution_promising, improvement,reward_rl] = rbf_ls_arm([population,objs_acq], database, paras.func_exp);
                paras.save_rl = [paras.save_rl, reward_rl];
                paras.rbf_model = [paras.rbf_model;reward_rl];
            elseif id == 4 % {GP, EI}
                [solution_promising, improvement,reward_ge] = gp_ei_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                paras.save_ge = [paras.save_ge, reward_ge];
                paras.gp_model = [paras.gp_model;reward_ge];
            elseif id == 5 % {PRS prescreening}
                [solution_promising, improvement,reward_pp] = prs_pre_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                paras.save_pp = [paras.save_pp, reward_pp];
                paras.prs_model = [paras.prs_model;reward_pp];
            elseif id == 6 % {PRS, local search}
                [solution_promising, improvement,reward_pl] = prs_ls_arm([population,objs_acq], database, paras.func_exp);
                paras.save_pl = [paras.save_pl, reward_pl];
                paras.prs_model = [paras.prs_model;reward_pl];
            elseif id == 7 % {KNN, L1-exploitation}
                [solution_promising, improvement,reward_ki] = knn_eoi_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                paras.save_ki = [paras.save_ki, reward_ki];
                paras.knn_model = [paras.knn_model;reward_ki];
            elseif id == 8 % {KNN, L1-exploration}
                [solution_promising, improvement,reward_ko] = knn_eor_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                paras.save_ko = [paras.save_ko, reward_ko];
                paras.knn_model = [paras.knn_model;reward_ko];
            end
        else
            % elect the high-level arm and low-level arm by TL-UCB
            u_model_value = zeros(4,1);
            for i = 1:4
                if i == 1
                    sum_reward = paras.rbf_model;
                    q_value_m = mean(paras.rbf_model); % The value of RBF
                elseif i == 2
                    sum_reward = paras.gp_model; 
                    q_value_m = mean(paras.gp_model);  % The value of GP
                elseif i == 3
                    sum_reward = paras.prs_model; 
                    q_value_m = mean(paras.prs_model); % The value of PRS
                elseif i == 4
                    sum_reward = paras.knn_model; 
                    q_value_m = mean(paras.knn_model); % The value of KNN
                end
                u_model_value(i) = tl_ucb(sum_reward, id, q_value_m);  % calculate UCB value of high-level arms
            end
            idx = find(u_model_value == max(u_model_value));
            if length(idx) >= 1 % if the max value of UCB exceed one arm, then select the arm at random
                seq = randperm(size(idx, 1));
                idx = idx(seq(1));
            end
            if idx == 1 % RBF is selected at the first level arm
                u_rbf_value = zeros(2,1);
                for i = 1:2
                    if i == 1
                        sum_reward = paras.save_rp;
                        q_value = mean(paras.save_rp);
                    elseif i == 2
                        sum_reward = paras.save_rl; 
                        q_value = mean(paras.save_rl);
                    end
                    u_rbf_value(i) = tl_ucb(sum_reward, id, q_value);  % calculate UCB values of low-level arms associated with the RBF
                end
                idx2 = find(u_rbf_value == max(u_rbf_value));
                if length(idx2) >= 1 % if the max value of UCB exceed one arm, then select the arm at random
                    seq = randperm(size(idx2, 1));
                    idx2 = idx2(seq(1));
                end
                if idx2 == 1 % {RBF, prescreening}
                    [solution_promising, improvement,reward_rp] = rbf_pre_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                    paras.save_rp = [paras.save_rp, reward_rp];
                    paras.rbf_model = [paras.rbf_model;reward_rp];
                elseif idx2 == 2 % {RBF, local search}
                    [solution_promising, improvement,reward_rl] = rbf_ls_arm([population,objs_acq], database, paras.func_exp);
                    paras.save_rl = [paras.save_rl, reward_rl];
                    paras.rbf_model = [paras.rbf_model;reward_rl];
                end
            elseif idx == 2 % GP is selected to be the first level arm
                u_gp_value = zeros(2,1);
                for i = 1:2
                    if i ==1
                        sum_reward = paras.save_gl;
                        q_value =  mean(paras.save_gl);
                    elseif i == 2
                        sum_reward = paras.save_ge;
                        q_value =  mean(paras.save_ge);
                    end
                    u_gp_value(i) = tl_ucb(sum_reward, id, q_value);  % calculate UCB values of low-level arms associated with the gp
                end
                idx2 = find(u_gp_value == max(u_gp_value));
                if length(idx2) >= 1
                    seq = randperm(size(idx2, 1));
                    idx2 = idx2(seq(1));
                end
                if idx2 == 1 % {GP, LCB}
                    [solution_promising, improvement,reward_gl] = gp_lcb_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                    paras.save_gl = [paras.save_gl, reward_gl];
                    paras.gp_model = [paras.gp_model;reward_gl];
                elseif idx2 == 2 % {GP, EI}
                    [solution_promising, improvement,reward_ge] = gp_ei_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                    paras.save_ge = [paras.save_ge, reward_ge];
                    paras.gp_model = [paras.gp_model;reward_ge];
                end
            elseif idx == 3 % PRS is selected to be the first level arm
                u_prs_value = zeros(2,1);
                for i = 1:2
                    if i ==1
                        sum_reward = paras.save_pp;
                        q_value =  mean(paras.save_pp);
                    elseif i == 2
                        sum_reward = paras.save_pl;
                        q_value =  mean(paras.save_pl);
                    end
                    u_prs_value(i) = tl_ucb(sum_reward, id, q_value);  % calculate UCB values of low-level arms associated with the prs
                end
                idx2 = find(u_prs_value == max(u_prs_value));
                if length(idx2) >= 1
                    seq = randperm(size(idx2, 1));
                    idx2 = idx2(seq(1));
                end
                if idx2 == 1 % {PRS, Prescreening}
                    [solution_promising, improvement,reward_pp] = prs_pre_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                    paras.save_pp = [paras.save_pp, reward_pp];
                    paras.prs_model = [paras.prs_model;reward_pp];
                elseif idx2 == 2 % {PRS, local search}
                    [solution_promising, improvement,reward_pl] = prs_ls_arm([population,objs_acq], database, paras.func_exp);
                    paras.save_pl = [paras.save_pl, reward_pl];
                    paras.prs_model = [paras.prs_model;reward_pl];
                end
            elseif idx == 4 % KNN is selected to be the first level arm
                u_knn_value = zeros(2,1);
                for i = 1:2
                    if i ==1
                        sum_reward = paras.save_ki;
                        q_value =  mean(paras.save_ki);
                    elseif i == 2
                        sum_reward = paras.save_ko;
                        q_value =  mean(paras.save_ko);
                    end
                    u_knn_value(i) = tl_ucb(sum_reward, id, q_value);  % calculate UCB values of low-level arms associated with the prs
                end
                idx2 = find(u_knn_value == max(u_knn_value));
                if length(idx2) >= 1
                    seq = randperm(size(idx2, 1));
                    idx2 = idx2(seq(1));
                end
                if idx2 == 1 % {KNN, L1-exploitation}
                    [solution_promising, improvement,reward_ki] = knn_eoi_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                    paras.save_ki = [paras.save_ki, reward_ki];
                    paras.knn_model = [paras.knn_model;reward_ki];
                elseif idx2 == 2 % {KNN, L1-exploration}
                    [solution_promising, improvement,reward_ko] = knn_eor_arm([population,objs_acq], population_offspring, database, paras.func_exp);
                    paras.save_ko = [paras.save_ko, reward_ko];
                    paras.knn_model = [paras.knn_model;reward_ko];
                end
            end
        end
        %********************Model and Infill Criterion Auto-Configuration********************%
end