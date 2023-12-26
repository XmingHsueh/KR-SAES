function [database,objs,transfer_states,paras]= sas_ckt(fun,lb,ub,FEsMax,name_saea,knowledge_base,paras)

%**********************************************************************%
%************************Phase 1: Initialization************************%
%**********************************************************************%
if isempty(name_saea)
    name_saea = 'BO-LCB';
end
num_sources = length(knowledge_base);
interval_transfer = paras.interval;
num_initial_solutions = size(paras.warm_up,1);
solutions_initial = paras.warm_up;
objs_initial = zeros(num_initial_solutions,1);
for i = 1:num_initial_solutions
    objs_initial(i) = fun(solutions_initial(i,:));
end
FEsUsed = num_initial_solutions;
database = [solutions_initial, objs_initial];
objs = sort(objs_initial,'descend');
dim = length(lb);
transfer_states = [];
if strcmp(name_saea,'AutoSAEA')
    load('optimizers\auto-saea\autosaea_paras.mat');
    options_saea.func_exp = fun;
else
    options_saea = [];
end

% build the surrogate models of k source tasks
surrogates_source = struct; 
n_count = 0;
while n_count < num_sources
    i = n_count + 1;
    try % the construction of GP may fail occasionally in a few releases
        surrogates_source(i).func = surrogate_model([knowledge_base(i).solutions,knowledge_base(i).objs],'gpr');
        n_count = n_count+1;
    catch
        continue;
    end
end

% execute the task adaptation if needed
if paras.ada_flag == 0
    ada_vectors = zeros(num_sources,dim);
else
    solutions_target = (solutions_initial-repmat(lb,num_initial_solutions,1))./...
            (repmat(ub,num_initial_solutions,1)-repmat(lb,num_initial_solutions,1));
    objs_target = objs_initial;
    ada_vectors = solution_adaptation(surrogates_source,solutions_target,objs_target);
end

%**********************************************************************%
%************************Phase 2: Optimization***********************%
%**********************************************************************%
while FEsUsed < FEsMax
    
    % execute one single acquisition with the designated SAS engine
    [solution_in,improvement_in,options_saea] = acquisition_single(database,lb,ub,name_saea,options_saea);
    
    % competitive knowledge transfer
    if mod(FEsUsed,interval_transfer)==0 
        [solution_ex,improvement_ex] = knowledge_competition(database,lb,ub,...
            knowledge_base,surrogates_source,ada_vectors);
        if improvement_in >= improvement_ex
            solution_candidate = solution_in;
            transfer_states = [transfer_states;0];
        else
            solution_candidate = solution_ex;
            transfer_states = [transfer_states;1];
        end
    else
        solution_candidate = solution_in;
    end
    
    % evaluate the candidate solution and update the database
    epsilon = 5e-3;
    no_trials  = 50;
    scales = linspace(0.1,1,no_trials);
    c = 1;
    while min(max(abs(repmat(solution_candidate,size(database,1),1)-...
            database(:,1:end-1)),[],2)) < epsilon
        perturbation = scales((mod(c,no_trials)==0)*no_trials+(mod(c,no_trials)~=0)*...
            mod(c,no_trials))*(ub-lb).*(rand(1,dim)-0.5*ones(1,dim));
        solution_candidate = solution_candidate + perturbation;
        solution_candidate(solution_candidate>ub(1)) = ub(1);
        solution_candidate(solution_candidate<lb(1)) = lb(1);
        c = c+1;
    end
    obj_candidate = fun(solution_candidate);
    database = [database;solution_candidate obj_candidate];
    objs = [objs;(obj_candidate<=objs(end))*obj_candidate+...
        (obj_candidate>objs(end))*objs(end)];
    FEsUsed = FEsUsed+1;
    if paras.ada_flag == 0
        fprintf([name_saea,'-CKT-FEsUsed: %d, obj_best: %.2f\n'],FEsUsed,objs(end));
    else
        fprintf([name_saea,'-CKT-ADA-FEsUsed: %d, obj_best: %.2f\n'],FEsUsed,objs(end));
    end
    
end