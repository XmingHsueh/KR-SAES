function [database,objs]= SAES(fun,lb,ub,FEsMax,paras)

% initialization
surrogate_name = paras.model;
search_engine = paras.optimizer;
popsize = paras.popsize;
query_max = paras.query;
selector = paras.selection;
acquisition = paras.acquisition;
parent_strategy = paras.parent_pop;
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

while FEsUsed < FEsMax
    
    % construct the surrogate model
    flag = 0;
    while flag == 0
        try
            [func_surrogate,model_type] = surrogate_model(database,surrogate_name);
            flag = 1;
        catch
            continue;
        end
    end
    
    % surrogate-assisted evolutionary search
    no_queries = 0;
    [population_parent,objs_parent] = parent_pop(database,lb,ub,parent_strategy,popsize,func_surrogate,acquisition); % form a parental population
    while no_queries < query_max
        population_offspring = offspring_pop(population_parent,objs_parent,lb,ub,search_engine); % generate query solutions using the parental population
        objs_offspring = offspring_acquisition(func_surrogate,population_offspring,acquisition); % acquisite objective values of the query solutions
        [population_parent,objs_parent] = selection_pop(population_parent,objs_parent,population_offspring,objs_offspring,selector); % form the next parental population
        no_queries = no_queries+1;
    end
    [~,idx] = min(objs_parent);
    solution_candidate = population_parent(idx,:);
    
    % evaluate the candidate solution and update the database
    epsilon = 5e-3;
    no_trials  = 50;
    scales = linspace(0.1,1,no_trials);
    c = 1;
    while min(max(abs(repmat(solution_candidate,size(database,1),1)-database(:,1:end-1)),[],2)) < epsilon
        perturbation = scales((mod(c,no_trials)==0)*no_trials+(mod(c,no_trials)~=0)*mod(c,no_trials))*(ub-lb).*(rand(1,dim)-0.5*ones(1,dim));
        solution_candidate = solution_candidate + perturbation;
        solution_candidate(solution_candidate>ub(1)) = ub(1);
        solution_candidate(solution_candidate<lb(1)) = lb(1);
        c = c+1;
    end
    
    obj_candidate = fun(solution_candidate);
    database = [database;solution_candidate obj_candidate];
    objs = [objs;(obj_candidate<=objs(end))*obj_candidate+(obj_candidate>objs(end))*objs(end)];
    FEsUsed = FEsUsed+1;
    fprintf('/SAES/ FEsUsed: %d, obj_best: %.2f\n',FEsUsed,objs(end));

end