%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Association selected floes (M. Rabatel IR 09-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of file_name must be: 'file_name.h5'
%
function selection_floes_matching(filename)

    %% Preliminary infos:
    groupname               = 'floe_shapes';
    datasetname             = 'floe_states';
    
    filename_selec_floes    = 'out_partial.h5';
    
    sel_floes               = h5info(filename_selec_floes);
    floes                   = h5info(filename);
    nb_sel_floes            = size(sel_floes.Groups(1).Datasets,1);
    nb_floes                = size(floes.Groups(1).Datasets,1);
    
    SF_Shapes = cell(nb_sel_floes,1); F_Shapes = cell(nb_floes,1);
    SF_Center = zeros(nb_sel_floes,2); F_Center = zeros(nb_floes,2);
    
    member_name = strcat('/',datasetname);
    A = h5read(filename_selec_floes,member_name); % read in transposed ??!!
    SF_Center(:,1) = A(1,:,1); SF_Center(:,2) = A(2,:,1);
    A = h5read(filename,member_name); % read in transposed ??!!
    F_Center(:,1) = A(1,:,1); F_Center(:,2) = A(2,:,1);
    
    for i=1:nb_sel_floes
        member_name = strcat('/',groupname,'/',int2str(i-1));
        A = h5read(filename_selec_floes,member_name); % read in transposed ??!!
        SF_Shapes{i} = A';
    end
    
    for i=1:nb_floes
        member_name = strcat('/',groupname,'/',int2str(i-1));
        A = h5read(filename,member_name); % read in transposed ??!!
        F_Shapes{i} = A';
    end

    %% Comparison:
    list_idx = zeros(nb_sel_floes,1);
    for i=1:nb_sel_floes
        if mod(i,20)==0; fprintf('%.4g%%\n',i/nb_sel_floes*100); end
        X = F_Center(:,1) - SF_Center(i,1);
        Y = F_Center(:,2) - SF_Center(i,2);
        ZZ = abs(X)+abs(Y);
        [a,~] = find(ZZ<1e-7);
        list_idx(i) = a-1;
    end
    
    %% saving
    save('list_selected_floes','list_idx');
    h5create('selected_floes.h5','/selected_floes_ids',nb_sel_floes,'DataType','int64');
    h5write('selected_floes.h5','/selected_floes_ids',list_idx);
end
