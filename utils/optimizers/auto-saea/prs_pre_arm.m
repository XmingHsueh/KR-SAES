function [solution_candidate, improvement, reward] = prs_pre_arm(data_train, pop_offspring, data_all, func_exp)
% {PRS, prescreening}
func_surrogate = surrogate_model(data_train,'prs');
objs_offspring = func_surrogate(pop_offspring);
[~,idx]=min(objs_offspring);
solution_candidate = pop_offspring(idx, :);
improvement = min(data_train(:,end))-min(objs_offspring);

epsilon = 5e-3;
if min(max(abs(repmat(solution_candidate,size(data_all,1),1)-data_all(:,1:end-1)),[],2)) < epsilon
    reward = 0;
else
    obj_candidate = func_exp(solution_candidate); % We do not count this evaluation for the sake of calculating the reward outside
    reward = low_level_r(data_train(:,end),obj_candidate);
end