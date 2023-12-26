function [solution_candidate, improvement, reward] = knn_eoi_arm(data_train, pop_offspring, data_all, func_exp)
% {KNN, L1-exploitation}

level = 5;
sidx = 1:size(data_train,1); %ghx and ghf have been sorted
labels_train = ceil(sidx*level/size(data_train,1));
parents_level1 = data_train(find(labels_train==1),1:end-1);
mdl =ClassificationKNN.fit(data_train(:,1:end-1),labels_train,'Distance','minkowski');% knn classfier
mdl.DistParameter = 2;

labels_offspring = predict(mdl,pop_offspring); %rank 1
pop_select = find(labels_offspring==min(labels_offspring));
dist = zeros(length(pop_select),size(parents_level1,1));
for i = 1: length(pop_select)
    for j = 1:size(parents_level1,1)
        dist(i,j) =  sqrt(sum((pop_offspring(pop_select(i), :)-parents_level1(j, :)).^2, 2));
    end
end
dist_max = max(dist,[],2);
[~,idx] = min(dist_max);
solution_candidate = pop_offspring(pop_select(idx(1)),:);
improvement = mean(diff(data_train(1:ceil(size(data_train,1)/level),end)));

epsilon = 5e-3;
if min(max(abs(repmat(solution_candidate,size(data_all,1),1)-data_all(:,1:end-1)),[],2)) < epsilon
    reward = 0;
else
    obj_candidate = func_exp(solution_candidate); % We do not count this evaluation for the sake of calculating the reward outside
    reward = low_level_r(data_train(:,end),obj_candidate);
end