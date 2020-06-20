function file_name = simulink_model(adjacency_matrix)
% Based on the provided adjacency matrix make a network of HR Neurons. 
    
    % Determine laplacian matrix, i.e. matrix representation of a graph
    degree_matrix = diag(sum(adjacency_matrix,2));
    sigma = 1;
    laplacian_matrix =sigma*(degree_matrix-adjacency_matrix);
    
    % open new simulink model and load system containing a HR neuron with
    % sampled input
    HR_neuron = 'Sampled_HR_100t';
    new_system('BlockDiagram')
    open_system('BlockDiagram')
    load_system(HR_neuron)
    
    % determine number of Neurons and create a list of strings for names of
    % subsystems
    nNeurons = size(adjacency_matrix,1);
    NeuronLabels = strings(1,size(adjacency_matrix,1));
    Subsystems = strings(3,1);
    AddStr = strings();
    
    for i = 1:nNeurons
        
        NeuronLabels(i) = strcat("Y",string(i));
        Subsystems(i) = strcat('BlockDiagram','/', NeuronLabels(i));
        AddStr = append(AddStr,'-');  % makes a list of n - signs for add block
        
    end
    
    % add a mux for combining all the outputs of the separate neurons
    add_block('simulink/Commonly Used Blocks/Mux', 'BlockDiagram/MuxOutputs')
    set_param('BlockDiagram/MuxOutputs', 'Inputs', string(nNeurons))
    
    for i = 1:nNeurons
        
        % create copies of the sampled HR neuron as subsystems with name Y1, Y2, etc.
        add_block('built-in/Subsystem', Subsystems(i))
        Simulink.BlockDiagram.copyContentsToSubsystem(HR_neuron, Subsystems(i))
        
        % set the initial conditions of the neurons in the subsytem to be
        % unique
        set_param(strcat(Subsystems(i),'/','IntegratorY'), 'InitialCondition', strcat('y0','(',string(i),')'))
        set_param(strcat(Subsystems(i),'/','IntegratorZ1'), 'InitialCondition', strcat('z10','(',string(i),')'))
        set_param(strcat(Subsystems(i),'/','IntegratorZ2'), 'InitialCondition', strcat('z20','(',string(i),')'))
        
        % add an 'Add' block with nNeurons inputs with all - sign and
        % connect 'Add' block to the corresponding subsystem 
        add_block('simulink/Math Operations/Add', strcat('BlockDiagram/Add',NeuronLabels(i)))
        set_param(strcat('BlockDiagram/Add',NeuronLabels(i)), 'Inputs',AddStr)
        add_line('BlockDiagram',strcat('Add',NeuronLabels(i),'/','1'),strcat(NeuronLabels(i),'/','1'),'autorouting','on')
        
        % connect neuron output to 'Mux' block
        add_line('BlockDiagram',strcat(NeuronLabels(i),'/','1'),strcat('MuxOutputs','/',string(i)),'autorouting','on')
        
        % create 'Gain' blocks with value depending on the laplacian matrix
        % and connect the 'Gain' blocks to the 'Add' Blocks.
        for j = 1:nNeurons
            
            add_block('simulink/Commonly Used Blocks/Gain', strcat('BlockDiagram/Gain',NeuronLabels(j),NeuronLabels(i)))
            set_param(strcat('BlockDiagram/Gain',NeuronLabels(j),NeuronLabels(i)), 'Gain', string(laplacian_matrix(j,i)))
            add_line('BlockDiagram',strcat('Gain',NeuronLabels(j),NeuronLabels(i),'/','1'),strcat('Add',NeuronLabels(i),'/',string(j)))
            
        end

    end
    
    % create a 'To Workspace' block and connect it to the 'Mux' block which
    % combined all neuron outpus.
    add_block('simulink/Sinks/To Workspace', 'BlockDiagram/To WorkspaceY')
    set_param('BlockDiagram/To WorkspaceY','VariableName', 'y')
    set_param('BlockDiagram/To WorkspaceY','SaveFormat','Array')
    add_line('BlockDiagram','MuxOutputs/1','To WorkspaceY/1')
    
    % Connect neuron outputs to the corresponding 'Gain' blocks, GainY1Y2
    % is output from Y1 for Y2
    for i = 1:nNeurons
        
        for j = 1:nNeurons
            
            add_line('BlockDiagram',strcat(NeuronLabels(i),'/','1'),strcat('Gain',NeuronLabels(i),NeuronLabels(j),'/','1'))
            
        end
        
    end
    
    set_param('BlockDiagram','SimulationMode','Accelerator');
    
    % make a new folder for the simulink files
    
    if ~exist('Simulink Models', 'dir')
       mkdir('Simulink Models')
    end
    
    % create a file name, which indicates the amount of neurons. 
    file_name = strcat('Simulink Models\Synchronization_',string(nNeurons),'_HR_3D_Neurons_100t','.slx');
    
    % make sure the file name is uniqe
    i = 0;
    while isfile(file_name)
        
        file_name = strcat('Simulink Models\Synchronization_',string(nNeurons),'_HR_3D_Neurons_100t','_',string(i),'.slx');
        i = i + 1;
        
    end
        
    % Arrange the blocks and close the system saving it under the created
    % file name.
    Simulink.BlockDiagram.arrangeSystem('BlockDiagram')
    close_system('BlockDiagram',strcat(file_name))
end