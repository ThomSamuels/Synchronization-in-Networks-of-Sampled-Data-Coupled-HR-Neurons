clear; clc; close all;
tic

% For a given adjacency matrix, which represents a network of HR neurons
% coupled via sampled-data couplings, the algorithm will compute the full
% (and partial) synchronization regions.

% Many network topologies of small networks have been simulated, and only a
% few larger than 5 nodes, with none larger than 8 nodes. The algorithm
% worked fine for these networks, however many topologies have not yet been
% tested and, thus, there might be some bugs that are still undiscovered.

%% Input
% adjacency_matrix of the network for which the full synchronization area
% will be determined

adjacency_matrix = [0 1
                    1 0];
                
% 1 (Full Synchronization) or 2 (Partial Synchronization)
synchronization = 1;

% Enabling this will visualize the construction of the synchronization
% area, disabling this will just show you the final result at the end

intermediate_plots = true;

% If a simulink model of the network already exist, give the file path if
% none exist, leave empty string, "".

global Simulink_file

Simulink_file = "";
%% Parameters

% Simulation
% global parameters needed for simulation

global parameterStruct simulation_times
global z10 z20 y0 nNeurons
global sigma h manifolds partition

E = 3.25;                   % Parameter value for which Hindmarsh-Rose neuron operates in a chaotic bursting mode
parameterStruct.ReturnWorkspaceOutputs = 'on';
parameterStruct.Solver = 'ode23tb';             % Simulation solver
simulation_times = ["20","50","100"];           % Simulation times

% Synchronization
% global parameters needed for synchronization function

global n allowed_error synchronization_time

n = 5;                      % #set of intial conditions for a combination of sigma-h to determine synchronization
allowed_error = 1e-4;       % Allowed difference in y values between to neurons to be considered synchronous
synchronization_time = 1.5; % Time for which the error should be small for the neurons to be considered synchronous (approx twice the max time between bursts)

% Bisection method
% global parameters needed for bisection functions
global Tolerance_horizontal Tolerance_vertical
Tolerance_horizontal = 0.01;    % Minimal resolution 'horizontal' bisection method (fixing h)
Tolerance_vertical = 0.0001;    % Minimal resolution 'vertical' bisection method (fixing sigma)

% Synchronization area
m = 100;                % Number of validation points
A = 0.2;                % Maximum change in sigma value
B = 0.5;                % Increase of sigma interval if there is no synchronization at h = 0

%% Set-up

% create simulink model of the neuronal network
if Simulink_file == ""
    Simulink_file = simulink_model(adjacency_matrix);
end

% Laplacian matrix
laplacian_matrix = diag(sum(adjacency_matrix,2)) - adjacency_matrix;

% set of initial conditions
nNeurons = size(adjacency_matrix,1);
z10 = zeros(1,nNeurons);
z20 = zeros(1,nNeurons);
y0 = zeros(1,nNeurons);

% Determine which manifolds need to be considered.
if synchronization == 1
    manifolds = zeros(1,nNeurons);
else
    % compute all partial synchronization manifolds
    manifolds_total = compute_manifolds(adjacency_matrix,2);

    % find manifolds which have unique conditions that do not coincide with
    % those for synchronization.
    eigval_man = zeros(size(manifolds_total,1),2);
    L = laplacian_matrix;
    
    for j = 1:size(manifolds_total,1)
        
        count = zeros(1,length(unique(manifolds_total(j,:))));
        
        for i = 1:size(manifolds_total,2)
            count(manifolds_total(j,i)+1) = count(manifolds_total(j,i)+1) + 1;
        end
        
        ER = eye(size(manifolds_total,2));
        % Define Error matrix
        for i = 1:length(count)
            if count(i) > 1
                cluster = find(manifolds_total(j,:) == i - 1);
                for k = 2:length(cluster)
                    ER(cluster(k),cluster(1)) = 1;
                    ER(cluster(k),cluster(k)) = -1;
                end
            end
        end
        
        % Determine eigenvalues of ELE^-1
        P = ER*L*(ER\eye(size(manifolds_total,2)));   
        [row, col] = find(ER == -1);
        eigenval = round(eig(P(row,col)),3);
        eigval_man(j,:) = [min(eigenval) max(eigenval)];
    end
    
    % find list of manifolds that have unique conditions that do not
    % coincide with those for synchronization
    uniques = unique(eigval_man, 'rows');
    id_list = zeros(size(uniques,1),1);
    
    for i = 1:size(uniques,1)
        id_list(i) = find(eigval_man(:,1) == uniques(i,1) & eigval_man(:,2) == uniques(i,2), 1);
    end
    
    [id_list, I] = sort(id_list);
    manifolds = manifolds_total(id_list,:);
    uniques = uniques(I,:);
    list = unique(uniques(:,1));
end

% Empty lists to store boundary points of the synchronization regions.
X = zeros(size(manifolds,1),100)*NaN; % empty lists to store the computed boundary points
Y = zeros(size(manifolds,1),100)*NaN;
X_synch = zeros(1,m)*NaN;   % empty lists to store the validation sigma values for which there is synchronization
Y_synch = zeros(1,m)*NaN;   % empty lists to store the validation h values for which there is synchronization
X_nosynch = zeros(1,m)*NaN; % empty lists to store the validation sigma values for which there is no synchronization
Y_nosynch = zeros(1,m)*NaN; % empty lists to store the validation h values for which there is no synchronization

max_distance = 0.75;        % minimal distance between two consecutive points
min_sigma_diff = 0.05;      % minimal sigma difference between two points for additional measurement point

for partition = 1:size(manifolds,1)
    tic
    i = 1;                      % index to store the values
    %% Find sigma threshold value of the full synchronization area
    
    elapsed_time = toc;
    fprintf('%6.0fs: Determine sigma threshold value of the full synchronization area\n', elapsed_time)
        
    if partition == 1
        flag = true;
    elseif partition > 1
        if uniques(partition,1) > max(uniques(1:partition-1,1))
            flag = true;
        else
            flag = false;
        end
    end
    
    if flag
        
        h = 0;          % h value for finding sigma
        sigma = 0;
        synch = false;  % for while loop

        % find a sigma value for which there is full synchronization
        while ~synch

            % random sigma value, sigma interval increases stepwise
            sigma = sigma + B*rand();

            % determine synchronization for (sigma,h)
            synch = synchronization_function(adjacency_matrix,manifolds,partition,n);

        end

        % horizontal bisection method to find  sigma threshold value of the full synchronization area
        [a,b] = bisection_function_horizontal(0,sigma);

        % sigma threshold, minimal sigma value for full synchronization
        sigma_threshold = (a+b)/2;

        if partition == 1 
            sigma_threshold_full = b;
        end
        
        X(partition,i) = sigma_threshold; % points stored to plot
        Y(partition,i) = h;
        
        % plot computed boundary point
        if intermediate_plots
            figure(1)
            plot(X(partition,i),Y(partition,i), 'k.')
            hold on
        end
        i = i + 1;
    end

    %% Random increase of sigma and determine upper bound of synchronization area (sigma * h < gamma*)
    % first coarse outline of the full synchronization area

    elapsed_time = toc;
    fprintf('%6.0fs: Determine Coarse outline of the full synchronization area\n', elapsed_time)
    
    if flag
        sigma = b;    % start value of sigma
    else
        sigma = 0;
    end
    
    synch = false;
    
    % move along an existing boundary layer of the x-axis until there is
    % synchronition.
    while sigma < sigma_threshold_full + 3
        flag = true;      % for while loop
        if isnan(X(partition,1))
            while ~synch

                sigma = sigma + A*rand();
                
                % use existing boundary as lower bound
                for j = 1:partition-1
                    XX = X(j,~isnan(X(j,:)));
                    YY = Y(j,~isnan(Y(j,:)));
                    h_inter = interp1(XX,YY,sigma,'linear','extrap');
                    
                    if j == 1 && h_inter >= 0
                        h = h_inter;
                    elseif h_inter > h
                        h = h_inter;
                    elseif h_inter < 0 && h == 0
                        h = 0;
                    end
                end

                h = h + 2e-3*rand();

                synch = synchronization_function(adjacency_matrix,manifolds,partition,n);
                
                if sigma > sigma_threshold_full + 3
                    X(partition,i) = 0; 
                    flag = false;
                    break
                end

            end
        end
        
        % In case of synchronization, determine when there is no
        % synchronization anymore and use the bisection method to determine
        % the upper bound.
        while flag

            if i ~= 1
                sigma = sigma + A*rand();   % small increase in sigma
            end

            synch = true;               % for while loop
            h_lowerbound = 0;
            
            % use existing boundary as lower bound for the bisection method
            if partition == 1 || sigma < min(min(X(1:partition-1,:)))
                if (h - 0.25e-2) > 0 
                    h = h - 0.25e-2;
                else
                    h = 0;
                end
            else
                for j = 1:partition-1
                    XX = X(j,~isnan(X(j,:)));
                    YY = Y(j,~isnan(Y(j,:)));
                    h_inter = interp1(XX,YY,sigma);
                    if h_inter > h_lowerbound
                        h_lowerbound = h_inter;
                    end
                end
                if (h - 0.25e-2) > h_lowerbound
                    h = h - 0.25e-2;
                else
                    h = h_lowerbound;
                end
            end

            % find h value for which there is no synchronization
            while synch
                
                h = h + 2e-3*rand();

                % determine synchronization for (sigma,h)
                synch = synchronization_function(adjacency_matrix,manifolds,partition,1);

            end

            % vertical bisection method to find  upper bound of synchronization area
            [a,b] = bisection_function_vertical(h_lowerbound,h);
            
            if a - h_lowerbound < Tolerance_vertical
                break
            end

            X(partition,i) = sigma;       % store point for plotting
            Y(partition,i) = a;
            
            % termination condition, if h value has reduced again.
            if sigma > sigma_threshold_full + 3
                flag = false;
            else
                flag = true;
            end
            
            % plot the found boundary point
            if intermediate_plots
                figure(1)
                plot(X(partition,i),Y(partition,i), 'k.')
                hold on
            end
            
            i = i + 1;
        end
    end

    %% Increase grid resolution by adding additional measurement points between consecutive points which are far apart

    elapsed_time = toc;
    fprintf('%6.0fs: Increasing outline resolution of the full synchronization area\n', elapsed_time)

    start = 1;              % start at the first point
    flag = true;            % for while loop

    % if distance between two consecutive points is 
    while flag

        distance = zeros(1,sum(~isnan(X(partition,:)))-1);

        % distance between two cosecutive points
        for j = 1:length(distance)
            distance(j) = ((X(partition,j)-X(partition,j+1))^2+(100*(Y(partition,j)-Y(partition,j+1)))^2)^(1/2);
        end

        % if distance it too big, have an additinal measurement point at sigma exactly in between those points
        % once every point is considered, start at first point again until distance requirements are satisfied
        for j = start:length(distance)

            if distance(j) > max_distance && (X(partition,j+1)-X(partition,j)) > min_sigma_diff

                % sigma value of new point, exactly in between the two points
                sigma = (X(partition,j)+X(partition,j+1))/2;
                
                if min([Y(partition,j),Y(partition,j+1)]) - 0.25e-2 > 0  
                    h = min([Y(partition,j),Y(partition,j+1)]) - 0.25e-2;
                else 
                    h = 0;
                end

                synch = true; % for while loop
                p = 1;
                B = 1;
                
                % use existing boundary as lower bound
                if partition == 1
                        h_lowbound = 0;
                else
                    for k = 1:partition-1
                        XX = X(k,~isnan(X(1,:)));
                        YY = Y(k,~isnan(Y(1,:)));
                        h_inter = interp1(XX,YY,sigma);
                        if h_inter > h
                            h_lowbound = h_inter;
                            h = h_lowbound;
                        end
                    end
                end
                
                % find h value for which there is no synchronization
                while synch

                    h = h + 2e-3*rand();  

                    % check for synchronization
                    synch = synchronization_function(adjacency_matrix,manifolds,partition,1);
                end

                [a,b] = bisection_function_vertical(h_lowbound,h);

                % insert new points in list of points, sorted on increasing sigma
                X_old_nan = X(partition,:);
                Y_old_nan = Y(partition,:);
                X_old = X_old_nan(~isnan(X_old_nan));
                Y_old = Y_old_nan(~isnan(Y_old_nan));
                
                % insert new found point, such that the sigma values are
                % ordered.
                if sum(isnan(X(partition,:))) == 0
                    X(partition,1:sum(~isnan(X(partition,:)))+1) = [X_old(1:j) sigma X_old(j+1:end)];
                    Y(partition,1:sum(~isnan(Y(partition,:)))+1) = [Y_old(1:j) a Y_old(j+1:end)];
                    X(X == 0) = NaN;
                    Y(Y == 0) = NaN;
                else
                    X(partition,1:sum(~isnan(X(partition,:)))+1) = [X_old(1:j) sigma X_old(j+1:end)];
                    Y(partition,1:sum(~isnan(Y(partition,:)))+1) = [Y_old(1:j) a Y_old(j+1:end)];
                end

                if intermediate_plots
                    plot(X(partition,j+1),Y(partition,j+1), 'k.')
                end

                start = j + 1;
                break
            end

            % whenever the last point has been reached, start over again
            if j == length(distance)
                start = 1;
            end

        end

        k = 1;  % for while loop

        % if all distances between consecutive distance is smaller than max distance or sigma value smaller than min_sigma_diff
        % then the grid resolution is sufficient.
        while flag

            if distance(k) < max_distance || (X(partition,k+1)-X(partition,k)) < min_sigma_diff

                if k == length(distance)
                    flag = false;
                end

            else 
                break
            end

            k = k + 1;

        end
    end
end

%% validation

% plot the boundary of the full synchronization area with a continuous line
% through all points and if only the synchronization region is computed
% plot m validation points to validate the computed region

if synchronization == 1
    elapsed_time = toc;
    fprintf('%6.0fs: Validating the full synchronization area\n', elapsed_time)

    ii = 1;
    jj = 1;
    % random validation points validation points
    for l = 1:m

        sigma = max(X)*rand();
        h = ((floor(max(100*Y))+1)*1e-2)*rand();

        % check synchronization
        synch = synchronization_function(adjacency_matrix,manifolds,partition,n);

        if synch
            X_synch(ii) = sigma;
            Y_synch(ii) = h;

            if intermediate_plots
                plot(X_synch(ii),Y_synch(ii), 'g.')
                hold on
            end

            ii = ii + 1;

        else
            X_nosynch(jj) = sigma;
            Y_nosynch(jj) = h;

            if intermediate_plots
                plot(X_nosynch(jj),Y_nosynch(jj), 'r.')
                hold on
            end

            jj = jj + 1;
        end

    end
end

% plot final points based sorted array

plot(X',Y',...
   'k.:','LineWidth',1.5,...
   'MarkerEdgeColor','k',...
   'MarkerFaceColor','k',...
   'MarkerSize',12);
hold on

% plot features
xlabel('\sigma');
ylabel('h [s]');
axis([0 max(X) 0 (floor(max(100*Y))+1)*1e-2]);
y_labels = 0:1:(floor(max(100*Y))+1);
yticks(y_labels/100);
yticklabels({y_labels});

% visualization of the validation points
plot(X_nosynch,Y_nosynch,'r.','MarkerSize',8)
plot(X_synch,Y_synch,'g.','MarkerSize',8)

ax = gca;
ax.FontSize = 14;

% plot all boundaries in a single plot.
figure(size(manifolds,2) + 1)
hold on
plot(X',Y','k.:')

%%
elapsed_time = toc;
fprintf('%6.0fs: Finished\n', elapsed_time)