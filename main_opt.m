% Author: Xiaoming Xue
% Email: xminghsueh@gmail.com
%
% ------------
% Description:
% ------------
% This file is the entry point for running KR-SAES.
%
% ------------
% Reference:
% ------------
% X. Xue, L. Feng, Y. Hu, et al. "Knowledge Race-Enhanced Surrogate-Assisted Evolutionary 
% Search of Computationally Expensive Problems", Submitted for Peer Review.

%% Initialization
clc,clear
warning off;

task_families = {'Sphere','Ellipsoid','Schwefel','Quartic','Ackley','Rastrigin','Griewank','Levy'}; % eight task families
transfer_scenarios = {'a','e'}; % intra-family and inter-family transfers
xis = [0 0.7 1]; % the parameter xi that controls the optimum coverage
similarity_distributions = {'c','u','i','d'}; % four representative similarity distributions
k = 20; % the number of previously-solved source tasks
folder_problems = '.\benchmarks';
specifications = [1 1 1 1 50 k; % STOP 1
    2 2 1 2	25 k; % STOP 2
    3 1 1 3	30 k; % STOP 3
    4 2 1 4	50 k; % STOP 4
    5 1 3 3	25 k; % STOP 5
    6 2 3 2	50 k; % STOP 6
    7 1 2 3	25 k; % STOP 7
    8 2 2 4	30 k; % STOP 8
    1 1 3 1	25 k; % STOP 9
    6 2 2 1	30 k; % STOP 10
    5 1 2 1	50 k; % STOP 11
    2 2 3 1	50 k]; % STOP 12
num_initial_solutions = 50;
paras.model = 'GPR'; % GPR, GLM, EnTr
paras.optimizer = 'EA';
paras.popsize = 50;
paras.query = 1; % Preselection: 1, Iteration: 20
paras.selection = 'RouletteWheel';
paras.acquisition = 'Plain';
paras.parent_pop = 'Elite';
paras.gen_gap = 20;
paras.FEsMaxAda = 0;
paras_ada = paras;
paras_ada.FEsMaxAda = 10000;
FEsMax = 500;
runs = 30; % the number of independent runs
problem_list = 1:size(specifications,1);
num_problems = length(problem_list); % the number of individual benchmark problems to be optimized
h=waitbar(0,'Starting'); % process monitor
runs_total = num_problems*runs;
count = 0*num_problems*runs;


%% KR-SAES
for i = 1:length(problem_list)
    n = problem_list(i);
    % import the black-box STO problem to be solved
    stop_tbo = STEOP('func_target',task_families{specifications(n,1)},'trans_sce',...
        transfer_scenarios{specifications(n,2)},'xi',xis(specifications(n,3)),...
        'sim_distribution',similarity_distributions{specifications(n,4)},'dim',...
        specifications(n,5),'k',specifications(n,6),'mode','opt','folder_stops',folder_problems);
    target_task = stop_tbo.target_problem;
    fun = target_task.fnc;
    lb = target_task.lb;
    ub = target_task.ub;
    results_opt = struct;
    paras_ada.ada_vectors = [];
    for r = 1:runs
        paras.warm_up = lhsdesign_modified(num_initial_solutions,lb,ub);
        paras_ada.warm_up = paras.warm_up;
        [~,objs_saea]= SAEA(fun,lb,ub,FEsMax,paras);
        [~,objs_krsaea,output_kr]= KRSAES(fun,lb,ub,FEsMax,stop_tbo.knowledge_base,paras);
        [~,objs_kradasaea,output_krada,paras_ada]= KRSAES(fun,lb,ub,FEsMax,stop_tbo.knowledge_base,paras_ada);
        results_opt(r).objs = [objs_saea,objs_krsaea,objs_kradasaea];
        results_opt(r).extras = [output_kr,output_krada];
        count = count+1;
        waitbar(count/runs_total,h,sprintf('In progress: %.2f (%d-th run of the %d-th prob is done!)%%',...
            count/runs_total*100,r,n));
    end
    % save the results
    save(['.\results\',task_families{specifications(n,1)},'-T',...
        transfer_scenarios{specifications(n,2)},'-xi',num2str(xis(specifications(n,3))),...
        '-S',similarity_distributions{specifications(n,4)},'-d',num2str(specifications(n,5)),...
        '-k',num2str(specifications(n,6)),'-results.mat'],'results_opt');
end
close(h);