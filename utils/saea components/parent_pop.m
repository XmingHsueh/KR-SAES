function [population_parent,objs_parent] = parent_pop(database,lb,ub,parent_strategy,popsize,func_surrogate,acquisition)
population = database(:,1:end-1);
objs = database(:,end);

switch (parent_strategy)
    case 'New'
        population_parent = population(end-popsize+1:end,:);
        objs_parent = objs(end-popsize+1:end);
    case 'Elite'
        [~,idx] = sort(objs);
        population_parent = population(idx(1:popsize),:);
        objs_parent = objs(idx(1:popsize));
    case 'Elite+Random'
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
