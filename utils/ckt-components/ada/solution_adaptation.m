function ada_vectors = solution_adaptation(surrogates_source,solutions_target,objs_target)
num_sources = length(surrogates_source);
[~,dim] = size(solutions_target);
ada_vectors = zeros(num_sources,dim);
no_local_opt = 10;
for i = 1:num_sources
    func_ada = @(x)obj_ada(x,surrogates_source(i).func,solutions_target,objs_target);
    lb_ada = -1*ones(1,dim);
    ub_ada = 1*ones(1,dim);
    % perform the multi-start search based on the interior-point algorithm
    options = optimoptions('fmincon','Algorithm','interior-point','Display','off');
    x_ini = lhsdesign_modified(no_local_opt,lb_ada,ub_ada); x_ini(1,:) = zeros(1,dim);
    x_opt_ada = zeros(no_local_opt,dim);
    y_opt_ada = zeros(no_local_opt,1);
    for j = 1:no_local_opt
        [x_opt_ada(j,:),y_opt_ada(j)] = fmincon(func_ada,x_ini(j,:),[],[],[],[],lb_ada,ub_ada,[],options);
    end
    [~,idx] = min(y_opt_ada);
    ada_vectors(i,:) = x_opt_ada(idx,:);
    fprintf('The %d-th source task has been adapted!\n',i);
end

