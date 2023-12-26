function population_offspring = pop_offspring(population_parent,objs_parent,lb,ub,search_engine)
[popsize, dim] = size(population_parent);
population = (population_parent-repmat(lb,popsize,1))./(repmat(ub,popsize,1)-repmat(lb,popsize,1));
population_offspring = zeros(popsize,dim);

switch(search_engine)
    case 'ea'
        mu = 15;     % index of Simulated Binary Crossover (tunable)
        mum = 15;    % index of polynomial mutation
        probswap = 0.5; % probability of variable swap
        indorder = randperm(popsize);
        for i = 1:popsize/2
            p1 = indorder(i); % population_parent 1
            p2 = indorder(i+(popsize/2)); % population_parent 2
            u = rand(1,dim);
            cf = zeros(1,dim);
            cf(u<=0.5)=(2*u(u<=0.5)).^(1/(mu+1));
            cf(u>0.5)=(2*(1-u(u>0.5))).^(-1/(mu+1));

            % crossover
            pp1 = population(p1,:);
            pp2 = population(p2,:);
            child1 = 0.5*((1+cf).*pp1 + (1-cf).*pp2);
            child1(child1<0) = 0;
            child1(child1>1) = 1;
            pp1 = population(p2,:);
            pp2 = population(p1,:);
            child2 = 0.5*((1+cf).*pp1 + (1-cf).*pp2);
            child2(child2<0) = 0;
            child2(child2>1) = 1;

            % mutation
            temp1 = child1;
            for j=1:dim
                if rand(1)<1/dim
                    u=rand(1);
                    if u <= 0.5
                        del=(2*u)^(1/(1+mum)) - 1;
                        temp1(j)=child1(j) + del*(child1(j));
                    else
                        del= 1 - (2*(1-u))^(1/(1+mum));
                        temp1(j)=child1(j) + del*(1-child1(j));
                    end
                end
            end

            child1 = temp1;
            child1(child1<0) = 0;
            child1(child1>1) = 1;
            temp2 = child2;
            for j=1:dim
                if rand(1)<1/dim
                    u=rand(1);
                    if u <= 0.5
                        del=(2*u)^(1/(1+mum)) - 1;
                        temp2(j)=child2(j) + del*(child2(j));
                    else
                        del= 1 - (2*(1-u))^(1/(1+mum));
                        temp2(j)=child2(j) + del*(1-child2(j));
                    end
                end
            end
            child2 = temp2;
            child2(child2<0) = 0;
            child2(child2>1) = 1;

            % variable swap (uniform X)
            swap_indicator = (rand(1,dim) >= probswap);
            temp = child2(swap_indicator);
            child2(swap_indicator) = child1(swap_indicator);
            child1(swap_indicator) = temp;

            population_offspring(i,:) = lb+child1.*(ub-lb);
            population_offspring(i+popsize/2,:) = lb+child2.*(ub-lb);
        end
    case 'de-b'
        CR = 0.8;
        F = 0.5;
        [~,idx] = min(objs_parent);
        individual_best = population(idx,:);

        for i = 1:popsize
            % mutation
            r1 = randi([1,popsize]);
            while(r1==i)
                r1=randi([1,popsize]);
            end
            r2 = randi([1,popsize]);
            while (r2==i) || (r2==r1)
                r2 = randi([1,popsize]);
            end
            v_mutation = individual_best+F*(population(r1,:)-population(r2,:));

            % crossover
            individual_offspring = zeros(1,dim);
            p = randi([1,dim]);
            for j = 1:dim
                cr = rand;
                if cr<=CR || j==p
                    individual_offspring(j) = v_mutation(j);
                else
                    individual_offspring(j) = population(i,j);
                end
            end

            % boundary check
            individual_offspring(individual_offspring<0) = 0;
            individual_offspring(individual_offspring>1) = 1;
            population_offspring(i,:) = lb+individual_offspring.*(ub-lb);
        end
    case 'de-r'
        CR = 0.8;
        F = 0.5;
        [~,idx] = min(objs_parent);
        individual_best = population(idx,:);

        for i = 1:popsize
            % mutation
            r1 = randi([1,popsize]);
            while(r1==i)
                r1=randi([1,popsize]);
            end
            r2 = randi([1,popsize]);
            while (r2==i) || (r2==r1)
                r2 = randi([1,popsize]);
            end
            r3 = randi([1,popsize]);
            while (r3==i) || (r3==r2) || (r3==r1)
                r3 = randi([1,popsize]);
            end
            v_mutation = population(1,:)+F*(population(r2,:)-population(r3,:));

            % crossover
            individual_offspring = zeros(1,dim);
            p = randi([1,dim]);
            for j = 1:dim
                cr = rand;
                if cr<=CR || j==p
                    individual_offspring(j) = v_mutation(j);
                else
                    individual_offspring(j) = population(i,j);
                end
            end

            % boundary check
            individual_offspring(individual_offspring<0) = 0;
            individual_offspring(individual_offspring>1) = 1;
            population_offspring(i,:) = lb+individual_offspring.*(ub-lb);
        end
    case 'pso'

end