function synch = synchronization_function(adjacency_matrix,manifolds,partition,n)
% function to determine whether the network of HR neurons coupled via 
% sampled-data couplings will synchronize for a combination of (sigma-h)

global h allowed_error synchronization_time simulation_times parameterStruct Simulink_file
global z10 z20 y0

% (sigma-h) is considered to be in the synchronization region if for
% n random sets of initial value the diff in y value between all neurons
% diminishes. The spread in initial value, away from a nominal trajectory, 
% and simulation time will increase stepwise.

% preparation
count = zeros(1,length(unique(manifolds(partition,:))));
for i = 1:size(manifolds,2)
    count(manifolds(partition,i)+1) = count(manifolds(partition,i)+1) + 1;
end    

% empty lists to story values
signals = zeros(length(count(count > 1)),2);
signals2 = zeros(length(count(count < 2)),2);
synch_list = ones(1,2)*NaN;

for i = 1:n
    
    % simulation times
    if h == 0 || n == 1
        simulation_times = "100";
    else
        simulation_times = ["20","50","100"];
    end
    
    for j = 1:length(simulation_times)
        
        % tolerances
        RelTol = 6;
        AbsTol = 6;
        
        % simulation parameters
        parameterStruct.Stoptime = simulation_times(j);
        parameterStruct.RelTol = strcat('1e-', num2str(RelTol));
        parameterStruct.AbsTol = strcat('1e-', num2str(AbsTol));
        
        % perturbation from nominal trajectory
        perturbation = 0.0001*str2double(simulation_times(j))^2;
        if perturbation > 1
            perturbation = 1;
        end
        
        % initial conditions for all neurons
        for k = 1:size(adjacency_matrix,1)
            z10(k) = (1-perturbation)*-0.2916 + perturbation*((0.9307+0.613)*rand()-0.613);
            z20(k) = (1-perturbation)*-2.7873 + perturbation*((-2.5642+3.038)*rand()-3.038);
            y0(k) = (1-perturbation)*-2.1389 + perturbation*((0.8+2.27)*rand()-2.27);
        end
        
        % simulation
        SimOut = sim(Simulink_file,parameterStruct);
        
        % simulation outpus
        t = SimOut.tout;
        y = SimOut.y;
        
        % construction of T matrix which determines which signals need to
        % be compared to each other. Consistst of two steps, first nodes in
        % the same cluster need to be compared to each other and then each
        % cluster need to be compared to each other.
        
        % Preparation of matrix T
        v = 1;
        w = 1;
        T_length = 0;
        signals_2b_considered = count > 1;
        signals_not_2b_considered = count < 2;
        for k = 1:length(count)
            if signals_2b_considered(k) == true
                signals(v,:) = [k-1 count(k)];
                T_length = T_length + count(k)-1;
                v = v + 1;
            elseif signals_not_2b_considered(k) == true
                signals2(w,:) = [k-1 count(k)];
                w = w + 1;
            end
        end
        
        T_length2 = size(signals,1)*size(signals2,1);
        
        if size(signals(:,1),1) > 1
            T_length2 = T_length2 + size(nchoosek(signals(:,1),2),1);
        end
        
        T = zeros(T_length+T_length2,size(manifolds,2));
        w = 1;
        
        % Nodes is the same cluster.
        for k = 1:size(signals,1)
            % determine maximum error for each time step
            % max of [y1-y2 y1-y3 y1-y4 ....]' depending on size of network
            idx = find(manifolds(partition,:) == signals(k,1));
            for v = 1:length(idx)-1
                T(w,idx(1)) = 1;
                T(w,idx(v+1)) = -1;
                w = w + 1;
            end
        end
        
        error_list = ones(size(y,1),2)*NaN;
        
        if T_length == 1
            error_list(:,1) = (T(1:T_length,:) * y')';
        else
            errors = (T(1:T_length,:) * y');
            error_list(:,1) = max(abs(errors))';
        end
        
        % Nodes in different clusters
        
        % if all nodes belong to a cluster, then no synch between clusters
        if ~isempty(signals2)
            for k = 1:size(signals,1)
                first = find(manifolds(partition,:) == signals(k,1),1,'first');
                for l = 1:size(signals2,1)
                    second = find(manifolds(partition,:) == signals2(l,1),1, 'first');
                    T(w,first) = 1;
                    T(w,second) = -1;
                    w = w + 1;
                end
            end
        end
        
        % if there a groups with a single node in it, that one should not
        % synchronize with the other clusters
        if size(signals,1) > 1
            nk = nchoosek(signals(:,1),2);
            for k = 1:size(nk,1)
               first = find(manifolds(partition,:) == nk(k,1),1,'first');
               second = find(manifolds(partition,:) == nk(k,2),1, 'first');
               T(w,first) = 1;
               T(w,second) = -1;
            end
        end
        
        if T_length2 ~= 0
            if T_length2 == 1
                error_list(:,2) = (T(T_length+1:end,:) * y')';
            else
                errors = (T(T_length+1:end,:) * y');
                error_list(:,2) = max(abs(errors))';
            end
        end
        
        % The error between the nodes in the same cluster should vanish,
        % while the error between nodes in other clusters should persists.
        % Both are tested.
        for k = [1 2]
            
            error = error_list(:,k);
        
            v = 1;
            w = 0;

        % determine how many times the error is within the allowed error
            while v < length(error)
                if abs(error(v)) <= allowed_error
                    w = w + 1;
                    while abs(error(v)) <= allowed_error && v < length(error)
                        v = v + 1;
                    end    
                end
                v = v + 1;
            end

            synchronization_time_list = zeros(1,w-1);
            v = 1;
            w = 1;    

            % determine the time of each period within the allowed error
            while v <= length(error)

                if abs(error(v)) <= allowed_error

                    i_synchronization1 = v;

                    while abs(error(v)) <= allowed_error && v < length(error)
                      v = v + 1;
                    end

                    i_synchronization2 = v;

                    synchronization_time_list(w) = t(i_synchronization2)- t(i_synchronization1);

                    w = w + 1;
                end

                v = v + 1;
            end

            % there is only synchronization if the error remains for
            % synchronization_time seconds within the allowed error
            synchronization = synchronization_time_list(synchronization_time_list > synchronization_time);

            % label the error profile as one of the following three:
            % synchronization, no synchronization or peaks. If no
            % synchronization then output is instantly synch = false, when
            % peaks then increase the tolerances until there is either
            % synchronization or no synchronization

            % synchronization: if the synchronization list constains only one
            % value which corresponds to the last encountered period for which
            % the error is within the allowed error

            % peaks: whenever there are multiple instances of synchronization
            % which are all temporarily interrupted by small sudden isolated
            % peaks or whenever there is a long period of synchronization at
            % the end which is not the only period of synchronization but
            % covers the majority of the simulation time

            % no synchronization: if the synchronization list is empty or
            % whenever is is neither synchronization nor peaks

            % no synchronization
            if isempty(synchronization)
                synch_list(k) = false;
                flag = true;
                break
            % synchronization
            elseif (synchronization(1) == synchronization_time_list(end))
                synch_list(k) = true;
                flag = false;
            % peaks
            elseif  ( ...
                    length(synchronization) < 4 && length(synchronization) >= 1 ...
                    )...
                    || ...
                    (...
                    ((synchronization(end) == max(synchronization)) && synchronization(end) == synchronization_time_list(end))...
                    )
                synch_list(k) = 2;
                flag = false;

                AbsTol = AbsTol + 1;
                RelTol = RelTol + 1;
            % no synchronization
            else
                synch_list(k) = false;
                flag = true;
                break
            end

        % increase relative and absolute tolerances whenever the error is
        % labeled as peaks. Simulate again for the same set of initial
        % conditions. Same procedure as above.
            while synch_list(k) == 2
                while AbsTol < 9 || RelTol < 9

                    parameterStruct.RelTol = strcat('1e-', num2str(RelTol));
                    parameterStruct.AbsTol = strcat('1e-', num2str(AbsTol));

                    % simulation
                    SimOut = sim(Simulink_file,parameterStruct);

                    % simulation outpus
                    t = SimOut.tout;
                    y = SimOut.y;

                    % determine maximum error for each time step
                    % max of [y1-y2 y1-y3 y1-y4 ....]' depending on size of network
                    
                    error_list = ones(size(y,1),2)*NaN;
                    
                    if T_length == 1
                        error_list(:,1) = (T(1:T_length,:) * y')';
                    else
                        errors = (T(1:T_length,:) * y');
                        error_list(:,1) = max(abs(errors))';
                    end

                    if T_length2 ~= 0
                        if T_length2 == 1
                            error_list(:,2) = (T(T_length+1:end,:) * y')';
                        else
                            errors = (T(T_length+1:end,:) * y');
                            error_list(:,2) = max(abs(errors))';
                        end
                    end
                    
                    error = error_list(:,k);

                    v = 1;
                    w = 0;

                    while v < length(error)
                        if abs(error(v)) <= allowed_error
                            w = w + 1;
                            while abs(error(v)) <= allowed_error && v < length(error)
                                v = v + 1;
                            end    
                        end
                        v = v + 1;
                    end

                    synchronization_time_list = zeros(1,w-1);
                    v = 1;
                    w = 1;    

                    while v <= length(error)

                        if abs(error(v)) <= allowed_error

                            i_synchronization1 = v;

                            while abs(error(v)) <= allowed_error && v < length(error)
                              v = v + 1;
                            end

                            i_synchronization2 = v;

                            synchronization_time_list(w) = t(i_synchronization2)- t(i_synchronization1);

                            w = w + 1;
                        end

                        v = v + 1;
                    end

                    synchronization = synchronization_time_list(synchronization_time_list > synchronization_time);

                    if isempty(synchronization)
                        synch_list(k) = false;
                        flag = true;
                        break
                    elseif (synchronization(1) == synchronization_time_list(end))
                        synch_list(k) = true;
                        flag = false;
                        break
                    % peaks
                    elseif ( ...
                            length(synchronization) < 4 && length(synchronization) >= 1 ...
                            )...
                            || ...
                            (...
                            ((synchronization(end) == max(synchronization)) && synchronization(end) == synchronization_time_list(end))...
                            )
                        if AbsTol < 8 && RelTol < 8
                            AbsTol = AbsTol + 1;
                            RelTol = RelTol + 1;
                            synch_list(k) = 2;
                            flag = false;
                        else
                            synch_list(k) = false;
                            flag = true;
                            break
                        end
                    else
                        synch_list(k) = false;
                        flag = true;
                        break
                    end
                end
            end
            if flag
                break
            end
        end
        
        if synch_list(1) == true && synch_list(2) == false
            synch = true;
        else
            synch = false;
            flag = true;
        end
        
        if synch == false
            break
        end
    end
    
    if flag
        break
    end
end