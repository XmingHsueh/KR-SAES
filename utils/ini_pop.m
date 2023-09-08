function [population_parent,objs_parent] = ini_pop(database,lb,ub,ini_strategy,popsize,func_surrogate)

population = database(:,1:end-1);
objs = database(:,end);

switch (ini_strategy)
    case 'new'
        population_parent = population(end-popsize+1:end,:);
        objs_parent = objs(end-popsize+1:end);
    case 'elite'
        [~,idx] = sort(objs);
        population_parent = population(idx(1:popsize),:);
        objs_parent = objs(idx(1:popsize));
    case 'elite&random'
        [~,idx] = sort(objs);
        population_parent = population(idx(1:floor(popsize/2)),:);
        objs_parent = objs(idx(1:floor(popsize/2)));
        population_r = lhsdesign_modified(popsize-floor(popsize/2),lb,ub);
        population_parent = [population_parent;population_r];
        objs_r = zeros(popsize-floor(popsize/2),1);
        for i = 1:popsize-floor(popsize/2)
            objs_r(i) = func_surrogate(population_r(i,:));
        end
        objs_parent = [objs_parent;objs_r];
end
