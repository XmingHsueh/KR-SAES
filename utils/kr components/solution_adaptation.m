function ada_vectors = solution_adaptation(surrogates_source,solutions_target,objs_target,FEsAdaMax)
num_sources = length(surrogates_source);
[~,dim] = size(solutions_target);
ada_vectors = zeros(num_sources,dim);
for i = 1:num_sources
    func_ada = @(x)obj_ada(x,surrogates_source(i).func,solutions_target,objs_target);
    lb_ada = -1*ones(1,dim);
    ub_ada = 1*ones(1,dim);
    paras.optimizer = 'EA';
    paras.popsize = 100;
    paras.selection = 'Truncation';
    paras.warm_up = lhsdesign_modified(paras.popsize,lb_ada,ub_ada);
    paras.warm_up(1,:) = zeros(1,dim);
    [~,~,ada_vectors(i,:)] = EA(func_ada,lb_ada,ub_ada,FEsAdaMax,paras);
    fprintf('The adaptation phase of KR-SAES (%d out of %d)\n',i,num_sources);
end
