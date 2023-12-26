function [solution_candidate, improvement, reward] = gp_ei_arm(data_train, pop_offspring, data_all, func_exp)
% {GP, EI}

try
    func_surrogate = surrogate_model(data_train,'gpr');
    obj_best = min(data_all(:,end));
    eis = zeros(size(pop_offspring,1),1);
    for i =1:size(pop_offspring,1)
        [obj_pre,obj_std] = func_surrogate(pop_offspring(i,:));
        eis(i) = (obj_best-obj_pre)*normcdf((obj_best-obj_pre)/obj_std)+obj_std*normpdf((obj_best-obj_pre)/obj_std);
    end
    [improvement,idx]=max(eis);
    solution_candidate = pop_offspring(idx, :);

    epsilon = 5e-3;
    if min(max(abs(repmat(solution_candidate,size(data_all,1),1)-data_all(:,1:end-1)),[],2)) < epsilon
        reward = 0;
    else
        obj_candidate = func_exp(solution_candidate); % We do not count this evaluation for the sake of calculating the reward outside
        reward = low_level_r(data_train(:,end),obj_candidate);
    end
catch
    reward = 0;
end