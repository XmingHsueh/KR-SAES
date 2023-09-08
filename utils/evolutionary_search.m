function [database,objs,x_best,y_best] = evolutionary_search(fun,lb,ub,FEsMax,paras)

popsize = paras.popsize;
search_engine = paras.optimizer;
selector = paras.selection;
solutions_initial = paras.warm_up;
num_initial_solutions = size(paras.warm_up,1);
objs_initial = zeros(num_initial_solutions,1);
for i = 1:num_initial_solutions
    objs_initial(i) = fun(solutions_initial(i,:));
end
FEsUsed = num_initial_solutions;
database = [solutions_initial, objs_initial];
objs = sort(objs_initial,'descend');
population_parent = solutions_initial;
objs_parent = objs_initial;

while FEsUsed < FEsMax

    population_offspring = offspring_pop(population_parent,objs_parent,lb,ub,search_engine);
    
    objs_offspring = zeros(popsize,1);
    for i = 1:popsize
        objs_offspring(i) = fun(population_offspring(i,:));
    end
    FEsUsed = FEsUsed+popsize;

    [population_parent,objs_parent] = selection_pop(population_parent,objs_parent,population_offspring,objs_offspring,selector);

    database = [database;population_offspring objs_offspring];
    objs = [objs;(min(objs_parent)<=objs(end))*min(objs_parent)+(min(objs_parent)>objs(end))*objs(end)];
    
end
[y_best,idx] = min(database(:,end));
x_best = database(idx,1:end-1);