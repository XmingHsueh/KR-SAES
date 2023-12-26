function [solution_candidate, improvement, reward] = prs_ls_arm(data_train, data_all, func_exp)
% {PRS, local search}
func_surrogate = surrogate_model(data_train,'prs');
% perform the evolutionary search on the surrogate model
lb_local = min(data_train(:,1:end-1));
ub_local = max(data_train(:,1:end-1));
query_max = 100;
popsize = 50;
selector = 'roulette_wheel';
acquisition = 'plain';
ini_strategy = 'random';
search_engine = 'ea';
no_queries = 0;
[population,objs_acq] = pop_ini(data_train,lb_local,ub_local,popsize,ini_strategy,func_surrogate,acquisition); % initialize the population
while no_queries < query_max
    population_offspring = pop_offspring(population,objs_acq,lb_local,ub_local,search_engine); % generate the offspring population
    objs_acq_offspring = pop_acquisition(func_surrogate,population_offspring,acquisition,[]); % acquisite the objective values of the offspring
    [population,objs_acq] = pop_selection(population,objs_acq,population_offspring,objs_acq_offspring,selector); % form the next parental population
    no_queries = no_queries+1;
end
[~,idx]=min(objs_acq);
solution_candidate = population(idx, :);
improvement = improvement_internal(func_surrogate,acquisition,data_train(:,1:end-1),solution_candidate,[]);

epsilon = 5e-3;
if min(max(abs(repmat(solution_candidate,size(data_all,1),1)-data_all(:,1:end-1)),[],2)) < epsilon
    reward = 0;
else
    obj_candidate = func_exp(solution_candidate); % We do not count this evaluation for the sake of calculating the reward outside
    reward = low_level_r(data_train(:,end),obj_candidate);
end