% Author: Xiaoming Xue
% Email: xminghsueh@gmail.com
%
% ------------
% Description:
% ------------
% Surrogate-Assisted Evolutionary Search with Competitive Knowledge Transfer (SEAS-CKT)
%
% ------------
% Input:
% ------------
% fun: the objective function of a target task
% lb: the lower bound of decision space
% ub: the upper bound of decision space
% FEsMax: the maximum function evaluations
% knowledge_base: a knowledge base with a number of previously-solved source tasks
% paras: other parameters
%
% ------------
% Output:
% ------------
% database: a database for storing the evaluated solutions on the target task
% objs: the best objective values found by SAESCKT
%
% ------------
% Reference:
% ------------
% X. Xue, L. Feng, Y. Hu, et al. "Surrogate-Assisted Evolutionary Search with Competitive
% Knowledge Transfer for Expensive Optimization", Submitted to Peer Review.

function [database,objs]= SAESCKT(fun,lb,ub,FEsMax,knowledge_base,paras)

% initialization
surrogate_name = paras.model; % surrogate model
search_engine = paras.optimizer; % evolutionary solver for optimizing the surrogate
popsize = paras.popsize; % population size
query_max = paras.query; % the number of iterated populations at each round of acquisition
selector = paras.selection; % envorimental selection
acquisition = paras.acquisition; % acquisition manner
ini_strategy = paras.ini_pop; % population initialization strategy
eva_delta = paras.gen_gap; % the evaluation interval for triggering the knowledge competition
FEsAda = paras.FEsMaxAda; % the number of surrogate-based evaluations for task adaptation
dim = length(lb); % dimensionality
num_sources = length(knowledge_base); % the number of source tasks

num_initial_solutions = size(paras.warm_up,1); % the number of solutions for initializing the target surrogate model
solutions_initial = paras.ini_lhs;
objs_initial = zeros(num_initial_solutions,1);
for i = 1:num_initial_solutions
    objs_initial(i) = fun(solutions_initial(i,:));
end
FEsUsed = num_initial_solutions;
database = [solutions_initial, objs_initial];
objs = sort(objs_initial,'descend');

% construct (or load) the source surrogate models
surrogates_source = struct;
n_count = 0;
while n_count < num_sources
    i = n_count + 1;
    try % avoid failed constructions
        surrogates_source(i).func = surrogate_model(knowledge_base(i).database,surrogate_name);
        n_count = n_count+1;
    catch
        continue;
    end
end

% perform the task adaptation if needed
if FEsAda == 0
    ada_vectors = zeros(num_sources,dim); % no adaptation for solutions
else
    solutions_target = (solutions_initial-repmat(lb,num_initial_solutions,1))./(repmat(ub,num_initial_solutions,1)-repmat(lb,num_initial_solutions,1));
    objs_target = objs_initial;
    ada_vectors = task_adaptation(surrogates_source,solutions_target,objs_target,FEsAda);
end

while FEsUsed < FEsMax

    % construct (or update) the target surrogate model
    flag = 0;
    while flag == 0
        try % avoid failed constructions
            func_surrogate = surrogate_model(database,surrogate_name);
            flag = 1;
        catch
            continue;
        end
    end

    % perform the surrogate-assisted evolutionary search on the target surrogate model
    no_queries = 0;
    [population_parent,objs_parent] = ini_pop(database,lb,ub,ini_strategy,popsize,func_surrogate); % form a parental population
    while no_queries < query_max
        population_offspring = offspring_pop(population_parent,objs_parent,lb,ub,search_engine); % generate some query solutions with the parental population
        objs_offspring = offspring_acquisition(func_surrogate,population_offspring,acquisition); % acquire objective values of the query solutions
        [population_parent,objs_parent] = selection_pop(population_parent,objs_parent,population_offspring,objs_offspring,selector); % form the next parental population
        no_queries = no_queries+1;
    end
    [~,idx] = min(objs_parent); % get the most promising solution found by SAES
    improvement_internal  = objs(end)-objs_parent(idx); % estimate the internal improvement
    solution_internal = population_parent(idx,:); % the internal solution

    % competitive knowledge transfer, i.e., knowledge competition
    if mod(FEsUsed,eva_delta)==0
        % get the source solution with its estimated external improvement for the target task
        [solution_external,improvement_external] = knowledge_competition(database,lb,ub,knowledge_base,surrogates_source,ada_vectors);
        if improvement_internal >= improvement_external
            solution_candidate = solution_internal;
        else
            solution_candidate = solution_external;
        end
    else
        solution_candidate = solution_internal;
    end

    % reject neighbouring solutions with a prespecified threshold
    epsilon = 5e-3; % the threshold
    z = 1;
    no_trials  = 50;
    scales = linspace(0.1,1,no_trials);
    while min(max(abs(repmat(solution_candidate,size(database,1),1)-database(:,1:end-1)),[],2)) < epsilon
        perturbation = scales((mod(z,no_trials)==0)*no_trials+(mod(z,no_trials)~=0)*mod(z,no_trials))*(ub-lb).*(rand(1,dim)-0.5*ones(1,dim));
        solution_candidate = solution_candidate + perturbation;
        solution_candidate(solution_candidate>ub) = ub(solution_candidate>ub);
        solution_candidate(solution_candidate<lb) = lb(solution_candidate<lb);
        z = z+1;
    end

    % evaluate the candidate solution and update the database
    obj_candidate = fun(solution_candidate);
    database = [database;solution_candidate obj_candidate];
    objs = [objs;(obj_candidate<=objs(end))*obj_candidate+(obj_candidate>objs(end))*objs(end)];
    FEsUsed = FEsUsed+1;
    fprintf('/SAES-CKT/ FEsUsed: %d, obj_best: %.2f\n',FEsUsed,objs(end));

end
