clc,clear
warning off;
task_families = {'Sphere','Ellipsoid','Schwefel','Quartic','Ackley','Rastrigin','Griewank','Levy'}; % eight task families
transfer_scenarios = {'a','e'}; % intra-family and inter-family transfers
similarity_distributions = {'h1','h2','m1','m2','m3','m4','l1','l2'}; % eight similarity distributions
k = 50; % the number of source tasks
specifications = [1 1 1 50 k; % STOP 1
    2 2 2 25 k; % STOP 2
    3 1 2 30 k; % STOP 3
    4 2 1 50 k; % STOP 4
    5 1 3 25 k; % STOP 5
    6 2 4 50 k; % STOP 6
    7 1 5 25 k; % STOP 7
    8 2 6 30 k; % STOP 8
    1 1 7 25 k; % STOP 9
    6 2 8 30 k; % STOP 10
    5 1 8 50 k; % STOP 11
    2 2 7 50 k]; % STOP 12
folder_problems = '.\benchmarks';
sas_names = {'BO-LCB','TLRBF','GL-SADE','DDEA-MESS','LSADE','AutoSAEA'}; % six backbone SAS optimizers
num_initial_solutions = 50; % the number of initial solutions
FEsMax = 500; % the maximum number of expensive evaluations
runs = 30; % the number of independent runs
problem_list = 1:size(specifications,1);
num_problems = length(problem_list); % the number of individual benchmark problems to be optimized
no_algorithms = length(sas_names);
h=waitbar(0,'Starting'); % process monitor
runs_total = no_algorithms*num_problems*runs;
count = 0*no_algorithms*num_problems*runs;
paras.interval = 20; % the transfer interval in terms of the true evaluation
paras.ada_flag = 0; % whether to use the task adaptation, 1: yes, 0: no

for a =1:no_algorithms
    for i = 1:length(problem_list)
        n = problem_list(i);
        % import a sequential transfer optimization problem to be solved
        stop_tbo = STOP('func_target',task_families{specifications(n,1)},...
            'trans_sce',transfer_scenarios{specifications(n,2)},...
            'sim_distribution',similarity_distributions{specifications(n,3)},...
            'dim',specifications(n,4),...
            'k',specifications(n,5),...
            'mode','opt',...
            'folder_stops',folder_problems);
        target_task = stop_tbo.target_problem;
        fun = target_task.fnc;
        lb = target_task.lb;
        ub = target_task.ub;
        results_opt = struct;
        for r = 1:runs
            paras.warm_up = lhsdesign_modified(num_initial_solutions,lb,ub);
            % [database_saea,objs_saea]= sas(fun,lb,ub,FEsMax,sas_names{a},paras);
            [database_ktsaea,objs_ktsaea,states_ktsaea]= sas_ckt(fun,lb,ub,FEsMax,...
                sas_names{a},stop_tbo.knowledge_base,paras);
            results_opt(r).objs = [objs_saea,objs_ktsaea];
            results_opt(r).states = states_ktsaea;
            count = count+1;
            waitbar(count/runs_total,h,sprintf('In progress: %.2f%% (alg %d'' %d-th run on prob %d is done!)',...
                count/runs_total*100,a,r,n));
        end
        % save the results
        save(['.\results\portability\',sas_names{a},'\',task_families{specifications(n,1)},'-T',transfer_scenarios{specifications(n,2)},...
            '-S',similarity_distributions{specifications(n,3)},'-d',num2str(specifications(n,4)),...
            '-k',num2str(specifications(n,5)),'-results.mat'],'results_opt');
    end
end
close(h);