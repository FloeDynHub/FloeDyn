%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computation of the floe distribution exposent (M. Rabatel IR 09-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of file_name must be: 'file_name.h5'
%
function calc_exp_size_distrib(filename)

    %% Preliminary infos:
    groupname               = 'floe_shapes';
        
    floes                   = h5info(filename);
    nb_floes                = size(floes.Groups(1).Datasets,1);
    
    F_Shapes = cell(nb_floes,1);
    
    for i=1:nb_floes
        member_name = strcat('/',groupname,'/',int2str(i-1));
        A = h5read(filename,member_name); % read en transposed ??!!
        F_Shapes{i} = A';
    end

    %% Computation of the polygon area:
    d = zeros(nb_floes,1);
    for i=1:nb_floes
        if mod(i,ceil(nb_floes/20))==0; fprintf('%.4g%%\n',i/nb_floes*100); end
        A = polyarea(F_Shapes{i}(:,1),F_Shapes{i}(:,2));
        % Deduction of the diameter of the equivalent circle:
        d(i) = 2*sqrt(A/pi);
    end
    
    %% plot of the distribution as a cumulative probability (>d):
    x = linspace(min(d)-min(d)/100,max(d)+max(d)/100,5000);
    P = zeros(size(x,2),1);
    for i=1:size(x,2)
        P(i) = size(find(d>x(i)),1);
    end
    loglog(x,P,'r')
    hold on
    grid on
    
    % comparison with powar law distribution:
    % ln(y) = -\alpha ln( x/s_min ) + b, with y the floe number from 1 to
    % nb_floes, x the floe size from s_min to s_max. So when y=nb_floes, 
    % x=s_max and so b = ln(nb_floes);
    xx = linspace(min(d),max(d),5000);
    alpha = 1.5;
    yy = exp(-alpha*log(xx./min(d))+log(nb_floes));
    
%     d_alpha = zeros(nb_floes,1);
%     
%     idx = 1:1:nb_floes;
%     for i=1:nb_floes
%         d_alpha(i) = exp( -1/alpha*log(i) + log(max(d)) );
%     end
%     x_a = linspace(min(d_alpha)-min(d_alpha)/100,max(d_alpha)+max(d_alpha)/100,5000);
%     P_a = zeros(size(x_a,2),1);
%     for i=1:size(x_a,2)
%         P_a(i) = size(find(d_alpha>x_a(i)),1);
%     end
%     loglog(idx,d_alpha,'r')
    loglog(xx,yy,'k')
    
%     xx = 100:1:1e4;
%     yy = exp(-alpha.*log(xx));
%     plot(xx,yy,'--r')
    
    %% saving
    save('list_selected_floes','list_idx');
    h5create('selected_floes.h5','/selected_floes_ids',nb_sel_floes,'DataType','int64');
    h5write('selected_floes.h5','/selected_floes_ids',list_idx);
end
