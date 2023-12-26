function imp = improvement_internal(func_surrogate,acquisition,x_database,x_promising,FEsUsed)

switch(acquisition)
    case 'plain'
        objs_acq_database = func_surrogate(x_database);
        objs_acq_promising = func_surrogate(x_promising);
    case 'lcb'
        w = 2;
        [objs_pre_database,std_database] = func_surrogate(x_database);
        objs_acq_database = objs_pre_database-w*std_database;
        [objs_pre_promising,std_promising] = func_surrogate(x_promising);
        objs_acq_promising = objs_pre_promising-w*std_promising;
    case 'lcb-decay'
        w = 2-2*(1/(1+exp(5-20*FEsUsed/500)));
        [objs_pre_database,std_database] = func_surrogate(x_database);
        objs_acq_database = objs_pre_database-w*std_database;
        [objs_pre_promising,std_promising] = func_surrogate(x_promising);
        objs_acq_promising = objs_pre_promising-w*std_promising;
end
imp = min(objs_acq_database)-objs_acq_promising;