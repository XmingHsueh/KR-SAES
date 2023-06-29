function objs_offspring = offspring_acquisition(func_surrogate,population_offspring,acquisition)
popsize = size(population_offspring,1);
objs_offspring = zeros(popsize,1);

switch(acquisition)
    case 'Plain'
        for i = 1:popsize
            objs_offspring(i) = func_surrogate(population_offspring(i,:));
        end
    case 'LCB'
        % null
end