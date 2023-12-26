function [population,objs_acq] = pop_ini(database,lb,ub,popsize,ini_strategy,func_surrogate,acquisition)
population = database(:,1:end-1);
objs = database(:,end);

switch (ini_strategy)
    case 'new'
        population = population(end-popsize+1:end,:);
        objs_acq = objs(end-popsize+1:end);
    case 'elite'
        [~,idx] = sort(objs);
        population = population(idx(1:popsize),:);
        objs_acq = objs(idx(1:popsize));
    case 'random'
        population = lhsdesign_modified(popsize,lb,ub);
        switch (acquisition)
            case 'plain'
                objs_acq = func_surrogate(population);
            case 'lcb'
                w = 2;
                [objs_pre,objs_std] = func_surrogate(population);
                objs_acq = objs_pre-w*objs_std;
            case 'lcb-decay'
                w = 2-2*(1/(1+exp(5-20*size(database,1)/500)));
                [objs_pre,objs_std] = func_surrogate(population);
                objs_acq = objs_pre-w*objs_std;
        end
    case 'elite+random'
        [~,idx] = sort(objs);
        population = [population(idx(1:floor(popsize/2)),:);...
            lhsdesign_modified(popsize-floor(popsize/2),lb,ub)];
        switch (acquisition)
            case 'plain'
                objs_acq = func_surrogate(population);
            case 'lcb'
                w = 2;
                [objs_pre,objs_std] = func_surrogate(population);
                objs_acq = objs_pre-w*objs_std;
            case 'lcb-decay'
                w = 2-2*(1/(1+exp(5-20*size(database,1)/500)));
                [objs_pre,objs_std] = func_surrogate(population);
                objs_acq = objs_pre-w*objs_std;
        end
    case 'best+randperm'
        [~,idx] = min(objs);
        seq = randperm(length(objs));
        population = [population(idx,:);population(seq(1:popsize-1),:)];
        objs_acq = [objs(idx);objs(seq(1:popsize-1))];
end