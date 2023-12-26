function [population,objs] = pop_selection(population_parent,objs_parent,population_offspring,objs_offspring,selector)
popsize = size(population_parent,1);

switch(selector)
    case 'truncation'
        objs_total = [objs_parent;objs_offspring];
        population_total = [population_parent;population_offspring];
        [~,idx] = sort(objs_total);
        population = population_total(idx(1:popsize),:);
        objs = objs_total(idx(1:popsize));
    case 'roulette_wheel'
        objs_total = [objs_parent;objs_offspring];
        population_total = [population_parent;population_offspring];
        fit_total = objs_total - min(min(objs_total),0) + 1e-6;
        fit_total = cumsum(1./fit_total);
        fit_total = fit_total./max(fit_total);
        idx = arrayfun(@(S)find(rand<=fit_total,1),1:popsize);
        population = population_total(idx,:);
        objs = objs_total(idx);
    case 'comparison'
        population = zeros(popsize,size(population_parent,2));
        objs = zeros(size(population_parent,1),1);
        for i = 1:size(population_parent,1)
            if objs_parent(i)<objs_offspring(i)
                population(i,:) = population_parent(i,:);
                objs(i) = objs_parent(i);
            else
                population(i,:) = population_offspring(i,:);
                objs(i) = objs_offspring(i);
            end
        end
end
