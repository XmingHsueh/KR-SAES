function [solution_candidate, improvement, reward] = gp_lcb_arm(data_train, pop_offspring, data_all, func_exp)
% {GP, LCB}

try
    func_surrogate = surrogate_model(data_train,'gpr');
    acquisition = 'lcb';
    objs_acq_offspring = pop_acquisition(func_surrogate,pop_offspring,acquisition,[]);
    [~,idx]=min(objs_acq_offspring);
    solution_candidate = pop_offspring(idx, :);
    improvement = improvement_internal(func_surrogate,acquisition,data_train(:,1:end-1),solution_candidate,[]);
    
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

