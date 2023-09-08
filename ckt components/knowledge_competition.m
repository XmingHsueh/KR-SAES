function [solution_externel,imp_max] = knowledge_competition(database,lb,ub,knowledge_base,surrogates_source,ada_vectors)

num_sources = length(knowledge_base);
dim = length(lb);
solutions_target = (database(:,1:end-1)-repmat(lb,size(database,1),1))./(repmat(ub,size(database,1),1)-repmat(lb,size(database,1),1));
objs_target = database(:,end);
solutions = zeros(num_sources,dim);
improvements = zeros(num_sources,1);
similarities = zeros(num_sources,1);
ranks_target = zeros(length(objs_target),1);

for i = 1:length(objs_target)
    ranks_target(i) = sum(objs_target<objs_target(i))+1;
end
for i = 1:num_sources
    solutions_source = knowledge_base(i).database(:,1:end-1);
    objs_source = knowledge_base(i).database(:,end);
    objs_val = zeros(length(objs_target),1);
    for j = 1:length(objs_target)
        objs_val(j) = surrogates_source(i).func(solutions_target(j,:)-ada_vectors(i,:));
    end
    ranks_val = zeros(length(objs_target),1);
    for j = 1:length(objs_target)
        ranks_val(j) = sum(objs_val<objs_val(j))+1;
    end
    regularization = 1-max(abs(ada_vectors(i,:)));
    similarities(i) = regularization*mySpearman(ranks_target,ranks_val);
    improvements(i) = improvement_estimation(sort(objs_source,'descend'),sort(objs_val,'descend'),sort(objs_target,'descend'),similarities(i));
    [~,idx] = min(objs_source);
    x_ada_source = (solutions_source(idx,:)+ada_vectors(i,:));
    x_ada_source(x_ada_source<0) = 0;
    x_ada_source(x_ada_source>1) = 1;
    solutions(i,:) = lb+x_ada_source.*(ub-lb);
end
[imp_max,idx] = max(improvements);
solution_externel = solutions(idx,:);