function [database,objs]= sas(fun,lb,ub,FEsMax,name_saea,paras)

%**********************************************************************%
%************************Phase 1: Initialization************************%
%**********************************************************************%
if isempty(name_saea)
    name_saea = 'BO-LCB';
end
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
if strcmp(name_saea,'AutoSAEA')
    load('optimizers\auto-saea\autosaea_paras.mat');
    options_saea.func_exp = fun;
else
    options_saea = [];
end

%**********************************************************************%
%************************Phase 2: Optimization***********************%
%**********************************************************************%
while FEsUsed < FEsMax
    
    % execute one single acquisition with the designated SAS
    [solution_candidate,~,options_saea] = acquisition_single(database,lb,ub,name_saea,options_saea);
    
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
    fprintf([name_saea,'-FEsUsed: %d, obj_best: %.2f\n'],FEsUsed,objs(end));
    
end