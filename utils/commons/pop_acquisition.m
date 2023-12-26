function objs_acq_offspring = pop_acquisition(func_surrogate,population_offspring,acquisition,FEsUsed)

switch(acquisition)
    case 'plain'
        objs_acq_offspring = func_surrogate(population_offspring);
    case 'lcb'
        w = 2;
        [objs_pre,objs_std] = func_surrogate(population_offspring);
        objs_acq_offspring = objs_pre-w*objs_std;
    case 'lcb-decay'
        w = 2-2*(1/(1+exp(5-20*FEsUsed/500)));
        [objs_pre,objs_std] = func_surrogate(population_offspring);
        objs_acq_offspring = objs_pre-w*objs_std;
end