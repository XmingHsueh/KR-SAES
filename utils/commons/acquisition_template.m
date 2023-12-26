function [population,objs_acq] = acquisition_template(func_surrogate,lb,ub,configs)

% initialization
database = configs.database;
popsize = configs.popsize;
query_max = configs.query_max;
ini_strategy = configs.ini_strategy;
acquisition = configs.acquisition;
search_engine = configs.search_engine;
selector = configs.selector;
FEsUsed = size(database,1);

no_queries = 0;
[population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,func_surrogate,acquisition); % initialize the population
while no_queries < query_max
    population_offspring = pop_offspring(population,objs_acq,lb,ub,search_engine); % generate the offspring population
    objs_acq_offspring = pop_acquisition(func_surrogate,population_offspring,acquisition,FEsUsed); % acquisite the objective values of the offspring
    [population,objs_acq] = pop_selection(population,objs_acq,population_offspring,objs_acq_offspring,selector); % form the next parental population
    no_queries = no_queries+1;
end