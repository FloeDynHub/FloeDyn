%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trajectory Analyzes (M. Rabatel IR 09-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of file_name must be: 'file_name.h5'
%
function traj_analyze(filename)

    %% Preliminary infos:
    suf                     = split(filename,'out_partial_');
    suffix                  = split(suf{2},'.h5');
    mat_save_name           = strcat('traj_analyse_',suffix{1});
    dataname                = 'floe_states';
    areaname                = 'floe_shapes';
    kinename                = 'Kinetic Energy';
    timename                = 'time';
        
    floes                   = h5info(filename);
    try
        nb_time                 = floes.Datasets(3).Dataspace.Size(3);
        nb_floes                = floes.Datasets(3).Dataspace.Size(2);
    catch
        try
            nb_time                 = floes.Datasets(2).Dataspace.Size(3);
            nb_floes                = floes.Datasets(2).Dataspace.Size(2);
        catch
            nb_time                 = floes.Datasets(4).Dataspace.Size(3);
            nb_floes                = floes.Datasets(4).Dataspace.Size(2);
        end
    end
    assert(nb_time==60481 && (nb_floes==2000 || nb_floes==1 || nb_floes==500));
    
    fprintf('Kinetic Energy Recovering...\n');
    member_name = strcat('/',kinename);
    try
        KE = h5read(filename,member_name);
    catch
    end
    fprintf('Time Recovering...\n');
    member_name = strcat('/',timename);
    T = h5read(filename,member_name);
    
    % checking time
    figure(50)
    plot(1:nb_time,T);
    grid on
    
    %  | pos x | pos y | theta | speed x | speed y | rot | impulse received
    
    
    %% First step: removing mean trajectory (mass_center):
%     Traj_MC = h5read(filename,'/mass_center'); % be careful this is the whole pack barycenter 
%     Traj_MC = Traj_MC(:,1:nb_time);
    % preliminaries:
    F_Sizes = zeros(nb_floes,1); 
    Area = zeros(nb_floes,1);
    fprintf('Floe size computation...\n');
    for i=0:nb_floes-1
        member_name = strcat('/',areaname,'/',int2str(i));
        A = h5read(filename,member_name);
        Area(i+1) = polyarea(A(1,:),A(2,:));
        F_Sizes(i+1) = 2*sqrt( Area(i+1) / pi);
    end
       
    fprintf('Floe States Recovering...\n');
    member_name = strcat('/',dataname);
    floe_traj_x = h5read(filename,member_name,[1 1 1],[1 nb_floes nb_time]); 
    floe_traj_y = h5read(filename,member_name,[2 1 1],[1 nb_floes nb_time]); 
    floe_traj_x = reshape(floe_traj_x,nb_floes,nb_time);
    floe_traj_y = reshape(floe_traj_y,nb_floes,nb_time);
    
    % mass center of 2000 central floes
    if nb_floes~=1
        Mass_tot = sum(Area*917); % thickness = 1m!!
        MF_x = 917*Area.*floe_traj_x;
        MF_y = 917*Area.*floe_traj_y;
        Mass_C = [sum(MF_x)/Mass_tot ; sum(MF_y)/Mass_tot];
    else
        Mass_C = zeros(2,nb_time);
    end
    
    % removing mean trajectory
    fluct_traj_x = floe_traj_x-Mass_C(1,:);
    fluct_traj_y = floe_traj_y-Mass_C(2,:);
    
    %% floe comparison (size as the diameter of the disk and distance between floe pair)
    
    % biggest and smallest equal floe repartition:
    if nb_floes~=1
        Area_s = sort(Area);
        threshold = 99;
        idx_min_s = Area<Area_s(end-threshold);
        idx_max_s = Area>=Area_s(end-threshold);
    end
    
    % repartition by distances with each other:
    % fprintf('Warning!! taking a long time!! (>10m)');  
    %     % median initial distance + median pair floe size
    % dist_0 = zeros(1/2*nb_floes*(nb_floes-1),1);
    % pair_floe_S = zeros(1/2*nb_floes*(nb_floes-1),1);
    % p=1;
    % for i=1:nb_floes-1
    %     j=i+1;
    %     dist_0(p:p+nb_floes-j) = sqrt( (floe_traj_x(j:end,1)-floe_traj_x(i,1)).^2 +...
    %             (floe_traj_y(j:end,1)-floe_traj_y(i,1)).^2 );
    %     pair_floe_S(p:p+nb_floes-j) = F_Sizes(j:end) + F_Sizes(i);    
    %     p = p+nb_floes-j+1;   
    % end
    
    % med_dist0 = median(dist_0);
    % med_pair_f = median(pair_floe_S);
    
    % idx_min_d = dist_0 <= med_dist0;
    % idx_max_d = dist_0 >= med_dist0;
    % idx_min_pS = pair_floe_S <= med_pair_f;
    % idx_max_pS = pair_floe_S >= med_pair_f;

    % disp_r = zeros(1,nb_time-1); 
    % disp_r_min_d = zeros(1,nb_time-1);    disp_r_max_d = disp_r;
    % disp_r_min_S = disp_r;              disp_r_max_S = disp_r; 
%     k=1;
%     for i=1:nb_floes-1
%         if (mod(i,50)==0), fprintf('iteration: %.2g%%\n', i/nb_floes*100); end
%         for j=i+1:nb_floes
% %             v_temp = [floe_traj_x(j,:)-floe_traj_x(i,:) ; ...
% %                 floe_traj_y(j,:)-floe_traj_y(i,:)];
% %             disp_tmp = ( sqrt(v_temp(1,:).*v_temp(1,:) + v_temp(2,:).*v_temp(2,:)) ...
% %                 - sqrt(v_temp(1,1).*v_temp(1,1) + v_temp(2,1).*v_temp(2,1)) ).^2;
%             
%             r_i = [floe_traj_x(i,2:end)-floe_traj_x(i,1:end-1) ; ...
%                 floe_traj_y(i,2:end)-floe_traj_y(i,1:end-1)];
%             r_j = [floe_traj_x(j,2:end)-floe_traj_x(j,1:end-1) ; ...
%                 floe_traj_y(j,2:end)-floe_traj_y(j,1:end-1)];
%             disp_tmp = sum((r_i - r_j).^2,1);
%             
%             disp_r = disp_r + disp_tmp;
%             
%             if (idx_min_d(k))
%                 disp_r_min_d = disp_r_min_d + disp_tmp;
%             end 
%             if (idx_max_d(k))
%                 disp_r_max_d = disp_r_max_d + disp_tmp;
%             end
%             if (idx_min_pS(k))
%                 disp_r_min_S = disp_r_min_S + disp_tmp;
%             end 
%             if (idx_max_pS(k))
%                 disp_r_max_S = disp_r_max_S + disp_tmp;
%             end
%             k=k+1;
%         end
%     end
%     disp_r          = 2/(nb_floes*(nb_floes-1)-2)*disp_r;
%     disp_r_min_d    = 1/(size(find(idx_min_d==1),1)-1)*disp_r_min_d;
%     disp_r_max_d    = 1/(size(find(idx_max_d==1),1)-1)*disp_r_max_d;
%     disp_r_min_S    = 1/(size(find(idx_min_pS==1),1)-1)*disp_r_min_S;
%     disp_r_max_S    = 1/(size(find(idx_max_pS==1),1)-1)*disp_r_max_S;
  
    % data for statistics:
    temporal_window = 259200; % 3 days
    simu_left = 604800-temporal_window;
    % origin time along the trajectory:
    t_o = floor(linspace(1,(temporal_window/10),500)); % divided by 10 (min time_step)
    % time step:
    delta_t = [10:10:1e2 150:50:1e3 1500:500:1e4 12500:5000:simu_left]/10;
    % floe sets:
    if nb_floes>500
        floeSet = unique( ceil( nb_floes * rand(1,500) ) );
        floeSet_min_S = 1:1:2000; floeSet_min_S = floeSet_min_S(idx_min_s);
        floeSet_min_S = floeSet_min_S(unique( ceil( size(floeSet_min_S,2) * rand(1,500) ) ) );
    else
        floeSet = 1:1:nb_floes;
    end
    % pair floe:
    pair_l = unique( ceil( nb_floes * rand(1,50) ) ); 
    pair_r = unique( ceil( nb_floes * rand(1,250) ) );
    pair_r = setdiff(pair_r,pair_l);
    pair_floes = zeros(size(pair_r,2)*size(pair_l,2),2);
    k=1;
    for i=1:size(pair_l,2)
        pair_floes(k:k+size(pair_r,2)-1,1) = repmat(pair_l(i),size(pair_r,2),1);
        pair_floes(k:k+size(pair_r,2)-1,2) = pair_r;
        k = k+size(pair_r,2);
    end
    
    %% Variance computation of the distance to the origin:
    Diffu = zeros(size(delta_t));
    Diffu_min_s = zeros(size(delta_t));
    Diffu_max_s = zeros(size(delta_t));
    
    if nb_floes~=1
    fprintf('Variance computation of the distance for diffusion...\n')
    
    for i=1:size(delta_t,2)
        real_supp = find(t_o+delta_t(i)<=60481); 
        assert(length(real_supp) == length(t_o)); % the support must be the same for coherent comparison!!
        D_x = fluct_traj_x(floeSet,t_o(real_supp)+delta_t(i)) - fluct_traj_x(floeSet,t_o(real_supp));
        D_y = fluct_traj_y(floeSet,t_o(real_supp)+delta_t(i)) - fluct_traj_y(floeSet,t_o(real_supp));
        
        Vd = D_x.^2 + D_y.^2;
        vtemp = Vd(:);
        % fprintf('size of the sample: %d \n',size(vtemp,1));
        Diffu(i) = mean(vtemp);
    end
    end
    
    if nb_floes>500
    for i=1:size(delta_t,2)
        real_supp = find(t_o+delta_t(i)<=60481);
        assert(length(real_supp) == length(t_o)); % the support must be the same for coherent comparison!!
        D_x = fluct_traj_x(idx_max_s,t_o(real_supp)+delta_t(i)) - fluct_traj_x(idx_max_s,t_o(real_supp));
        D_y = fluct_traj_y(idx_max_s,t_o(real_supp)+delta_t(i)) - fluct_traj_y(idx_max_s,t_o(real_supp));
        
        Vd = D_x.^2 + D_y.^2;
        vtemp = Vd(:);
        % fprintf('size of the sample: %d \n',size(vtemp,1));
        Diffu_max_s(i) = mean(vtemp);

        %%% min
        D_x = fluct_traj_x(floeSet_min_S,t_o(real_supp)+delta_t(i)) - fluct_traj_x(floeSet_min_S,t_o(real_supp));
        D_y = fluct_traj_y(floeSet_min_S,t_o(real_supp)+delta_t(i)) - fluct_traj_y(floeSet_min_S,t_o(real_supp));
        
        Vd = D_x.^2 + D_y.^2;
        vtemp = Vd(:);
        % fprintf('size of the sample: %d \n',size(vtemp,1));
        Diffu_min_s(i) = mean(vtemp);
    end
    end
    
    %% autocorrelation:
    fprintf('Autocorrelation computation...\n')
    
    %speed:
    speed_x = 1/10*(fluct_traj_x(:,2:end) - fluct_traj_x(:,1:end-1) );
    speed_y = 1/10*(fluct_traj_y(:,2:end) - fluct_traj_y(:,1:end-1) );
    speed = sqrt(speed_x.*speed_x+speed_y.*speed_y);
% for i=1:nb_time-1
%         tmp_x = 1/10*(fluct_traj_x(:,i+1)-fluct_traj_x(:,i));
%         tmp_y = 1/10*(fluct_traj_y(:,i+1)-fluct_traj_y(:,i));
%         speed(:,i) = sqrt(tmp_x.*tmp_x+tmp_y.*tmp_y);
%     end
   
    % autocorrelation fct:
    % Comment: the mean mu and the variance var of the random process 
    % (speed any floe, any start time) are time dependent. That means, for
    % start time betweeen origin (t=0) to t=67500s (about 19h), mu and var
    % are respectively equal to 3.1e-3 and 1.2812e-5, and for start time
    % between t=177500(2j) to the end (4j) mu and var are respectively
    % equal to 7.4e-4 and 7.1e-7.
    % For a well-defined autocorrelation function, one need for mu and var
    % are time independent!
    % 
    % temporal window:
    temporal_window_a = 345600; % 3 days
    simu_left_a = 604800-temporal_window_a;
    % origin time along the trajectory:
    t_o_a = floor(linspace(1,(temporal_window_a/10),500)); % divided by 10 (min time_step)
    % time step:
    Frq = 48; % 48 per day -> 30 min -> 1800s
    day_in_s = 86400;
    delta_t_a = [0:day_in_s/Frq:simu_left_a]/10;
     
    ACF = ones(1,length(delta_t_a)); 
    ACF_min_s = zeros(1,length(delta_t_a)+1);
    ACF_max_s = zeros(1,length(delta_t_a)+1);
        
    N = length(t_o_a)*length(floeSet);
    mopt = zeros(1,length(delta_t_a)); vopt = mopt;
    V_o = speed(floeSet,t_o_a);
    vtemp = V_o(:);
    mo = mean(vtemp);
    vo = var(vtemp);
    for i=1:size(delta_t_a,2)
        real_supp = find(t_o_a+delta_t_a(i)<=60481);
        assert(length(real_supp) == length(t_o_a)); % the support must be the same for coherent comparison!!
        V_o_p_tau = speed(floeSet,t_o_a(real_supp)+delta_t_a(i));
        
        vtemp = V_o_p_tau(:);
        mopt(i) = mean(vtemp);
        vopt(i) = var(vtemp);
%         VV = (V_o-mo).*(V_o_p_tau-mopt(i));
        VV = (V_o-mo).*(V_o_p_tau-mo);
        
%         ACF(i+1) = sum(VV(:))/( (N-1) * sqrt(vo) * sqrt( vopt(i) ) );
        ACF(i) = sum(VV(:))/( (N-1) * vo );
    end
    
    if nb_floes>500
    N_max = length(t_o_a)*length(find(idx_max_s==1));
    N_min = length(t_o_a)*length(floeSet_min_S);
    mopt_min = mopt; vopt_min = vopt;
    mopt_max = mopt; vopt_max = vopt;
    
    V_o_max = speed(idx_max_s,t_o_a);
    vtemp = V_o_max(:);
    mo_max = mean(vtemp);
    vo_max = var(vtemp);
    
    V_o_min = speed(floeSet_min_S,t_o_a);
    vtemp = V_o_min(:);
    mo_min = mean(vtemp);
    vo_min = var(vtemp);
    
    for i=1:size(delta_t_a,2)
        real_supp = find(t_o_a+delta_t_a(i)<=60481);
        assert(length(real_supp) == length(t_o_a)); % the support must be the same for coherent comparison!!
        V_opt_max = speed(idx_max_s,t_o_a(real_supp)+delta_t_a(i));
        
        vtemp = V_opt_max(:);
        mopt_max(i) = mean(vtemp);
        vopt_max(i) = var(vtemp);
%         VV = (V_o_max-mo_max).*(V_opt_max-mopt_max(i));
        VV = (V_o_max-mo_max).*(V_opt_max-mo_max);
        
%         ACF_max_s(i+1) = sum(VV(:))/( (N_max-1) * sqrt(vo_max) * sqrt( vopt_max(i) ) );
        ACF_max_s(i+1) = sum(VV(:))/( (N_max-1) * vo_max );

        %%% min
        V_opt_min = speed(floeSet_min_S,t_o_a(real_supp)+delta_t_a(i));
        
        vtemp = V_opt_min(:);
        mopt_min(i) = mean(vtemp);
        vopt_min(i) = var(vtemp);
%         VV = (V_o_min-mo_min).*(V_opt_min-mopt_min(i));
        VV = (V_o_min-mo_min).*(V_opt_min-mo_min);
        
%         ACF_min_s(i+1) = sum(VV(:))/( (N_min-1) * sqrt(vo_min) * sqrt( vopt_min(i) ) );
        ACF_min_s(i+1) = sum(VV(:))/( (N_min-1) * vo_min );
    end
    end
        
    %% Power Spectrum Density (PSD) using Fourier transform:
    % from the autocorrelation with the ice floes velocity instead of the
    % fluctuation velocity (without the mean movement) in order to keep the
    % inertial oscillations:
    speed_x = 1/10*(floe_traj_x(floeSet,2:end) - floe_traj_x(floeSet,1:end-1) );
    speed_y = 1/10*(floe_traj_y(floeSet,2:end) - floe_traj_y(floeSet,1:end-1) );
    speed = sqrt(speed_x.*speed_x+speed_y.*speed_y);
    
    ACF_tmp = ones(1,length(delta_t_a)); 
    N = length(t_o_a)*length(floeSet);
    V_o = speed(:,t_o_a);
    vtemp = V_o(:);
    mo = mean(vtemp);
    vo = var(vtemp);
    for i=1:size(delta_t_a,2)
        real_supp = find(t_o_a+delta_t_a(i)<=60481);
        assert(length(real_supp) == length(t_o_a)); % the support must be the same for coherent comparison!!
        V_o_p_tau = speed(:,t_o_a(real_supp)+delta_t_a(i));
        
%         VV = (V_o-mo).*(V_o_p_tau-mopt(i));
        VV = (V_o-mo).*(V_o_p_tau-mo);
        
%         ACF(i+1) = sum(VV(:))/( (N-1) * sqrt(vo) * sqrt( vopt(i) ) );
        ACF_tmp(i) = sum(VV(:))/( (N-1) * vo );
    end
    
    L = length(delta_t_a);
    X = fft(ACF_tmp);
    PSD=X.*conj(X)/(L*L);
    PSD = PSD(1:ceil(L/2+1));
    f = Frq*(1:ceil(L/2+1))/L;
    clear ACF_tmp
    
    %% distribution of tangential and normal accelerations
    % from [Leveque and Naso 2014]
    % velocity computation 
    speed_x = 1/10*(fluct_traj_x(floeSet,2:end) - fluct_traj_x(floeSet,1:end-1) );
    speed_y = 1/10*(fluct_traj_y(floeSet,2:end) - fluct_traj_y(floeSet,1:end-1) );
    
    q = [0.5 1 1.5 2 2.5 3];
    acc = zeros(size(delta_t,2),length(q));
    for i=1:size(delta_t,2)
        if (mod(i,5)==0), fprintf('%.2g%%\n', i/size(delta_t,2)*100); end
        real_supp = find(t_o+delta_t(i)<=60481);

        V_x = 1/delta_t(i)*( speed_x(:,t_o(real_supp)+delta_t(i)) - speed_x(:,t_o(real_supp)) );
        V_y = 1/delta_t(i)*( speed_y(:,t_o(real_supp)+delta_t(i)) - speed_y(:,t_o(real_supp)) );

        for j=1:length(q)
            V = abs(V_x).^q(j)+abs(V_y).^q(j);
            vtemp = V(:);
            % fprintf('size of the sample: %d \n',size(vtemp,1));
            acc(i,j) = mean(vtemp);
        end
    end        

    % tau_dis = [10 50 100 250:250:1000 2000:1000:5000 10000:5000:30000 40000 50000 75000 100000]/10;
    
    % tang_a = zeros(size(tau_dis,2),nb_time);
    % norm_a = zeros(size(tau_dis,2),nb_time);
    % std_t  = zeros(size(tau_dis,2),1); std_n = std_t;
    % % frame:
    % for j=1:size(tau_dis,2)
    %     fprintf('iteration: %d\n',j);
    %     for i=1:nb_time-1-tau_dis(j)
    %         frame = [ fluct_traj_x(:,i+tau_dis(j)) - fluct_traj_x(:,i) , ...
    %             fluct_traj_y(:,i+tau_dis(j)) - fluct_traj_y(:,i) ];
    %         Nframe = sqrt( frame(:,1).^2 + frame(:,2).^2 );
    %         idx_mvt = Nframe~=0;
            
    %         if any(idx_mvt~=0)
    %             norm_frame = frame(idx_mvt,:)./Nframe(idx_mvt);
                
    %             tmp_xi = (fluct_traj_x(idx_mvt,i+1)-fluct_traj_x(idx_mvt,i));
    %             tmp_yi = (fluct_traj_y(idx_mvt,i+1)-fluct_traj_y(idx_mvt,i));
    %             tmp_xj = (fluct_traj_x(idx_mvt,i+tau_dis(j)+1)-fluct_traj_x(idx_mvt,i+tau_dis(j)));
    %             tmp_yj = (fluct_traj_y(idx_mvt,i+tau_dis(j)+1)-fluct_traj_y(idx_mvt,i+tau_dis(j)));
                
    %             speed_i = 1/10*[tmp_xi , tmp_yi];
    %             speed_j = 1/10*[tmp_xj , tmp_yj];
                
    %             gamma = 1/(10*tau_dis(j))*(speed_j-speed_i); % acceleration

    %             % decomposition in frame
    %             % B1 = [O;e1;e2] | B2 = [O';e1';e2']; O', e1' and e2' in B1
    %             % orthonal direct bases
    %             % e1' = a*e1 + b*e2 | e2' = c*e1 + d*e2
    %             % and a^2+b^2 = 1 and c = -b, d = a;
    %             %
    %             % P transition matrix from B2 to B1:
    %             %                 e1'    e2'
    %             %  P =    		  a      -b     e1
    %             %                 b       a     e2
    %             %
    %             % ex: 	e1' in B1
    %             %		(a,
    %             %		b)	= P	(1,
    %             %				0)
    %             %				e1' in B2
    %             %
    %             % e1 = P^-1 e'1;
    %             %          
    %             % P^-1 = a    b
    %             %        -b   a
    %             %
    %             % P^-1 x gamma(1)  =  g_1*nf_1  + g_2*nf_2
    %             %        gamma(2)     -g_1*nf_2 + g_2*nf_1
    %             decomp = [gamma(:,1).*norm_frame(:,1)+gamma(:,2).*norm_frame(:,2) , ...
    %                 gamma(:,2).*norm_frame(:,1)-gamma(:,1).*norm_frame(:,2)];

    %             % tangential acceleration
    %             tang_a(j,i) = mean(decomp(:,1));
    %             norm_a(j,i) = mean(decomp(:,2));
    %         end
    %     end
    %     std_t(j) = std(tang_a(j,1:nb_time-tau_dis(j)-1));
    %     std_n(j) = std(norm_a(j,1:nb_time-tau_dis(j)-1));
    % end
    
    %% Floe pair dispersion
    fprintf('Variance computation of the distance for floe pair dispersion...\n')
    
    Disper = zeros(size(delta_t));
%     for i=1:size(delta_t,2)
%         real_supp = find(t_o+delta_t(i)<=60481);
% 
%         D_x_1 = fluct_traj_x(pair_floes(:,1),t_o(real_supp)+delta_t(i)) - fluct_traj_x(pair_floes(:,1),t_o(real_supp));
%         D_y_1 = fluct_traj_y(pair_floes(:,1),t_o(real_supp)+delta_t(i)) - fluct_traj_y(pair_floes(:,1),t_o(real_supp));
%         
%         D_x_2 = fluct_traj_x(pair_floes(:,2),t_o(real_supp)+delta_t(i)) - fluct_traj_x(pair_floes(:,2),t_o(real_supp));
%         D_y_2 = fluct_traj_y(pair_floes(:,2),t_o(real_supp)+delta_t(i)) - fluct_traj_y(pair_floes(:,2),t_o(real_supp));
% 
%         R_x = D_x_1 - D_x_2; %R = r(i) - r(j)
%         R_y = D_y_1 - D_y_2; %R = r(i) - r(j)
% 
%         Vd = R_x.^2 + R_y.^2;
%         vtemp = Vd(:);
%         % fprintf('size of the sample: %d \n',size(vtemp,1));
%         Disper(i) = mean(vtemp);
%     end
    
    %% plot
%     fprintf('Plot...\n')
%     % switch suffix{1}
%     %     case 'HE6ob'
%     %         Vr = 402386; Vs =  [ -6.42729 0.298088 ]; Vo = [ 998926 -46328.7 ];
%     %     case 'cOlFJ'
%     %         Vr = 218938; Vs =  [ -1.82772 -9.25212 ]; Vo = [ 193801 981041 ];
%     %     case 'bSNvI'
%     %         Vr = 342577; Vs =  [ 3.73322 9.1694 ]; Vo = [ -377083 -926180 ];
%     %     case 'ik9dr' %='bSNvI'
%     %         Vr = 342577; Vs =  [ 3.73322 9.1694 ]; Vo = [ -377083 -926180 ];
%     %     case '9VTBp'
%     %         Vr = 374848; Vs =  [ -9.10191 0.813075 ]; Vo = [ 996034  -88975.9 ];
%     %     case '3sbVt' %='9VTBp'
%     %         Vr = 374848; Vs =  [ -9.10191 0.813075 ]; Vo = [ 996034  -88975.9 ]; 
%     %     case 'Ms3Dl'
%     %         Vr = 325639; Vs =  [ 1.26598 8.87528 ]; Vo = [ -141212  -989979 ];
%     %     case 'UkdVH' %='Ms3Dl'
%     %         Vr = 325639; Vs =  [ 1.26598 8.87528 ]; Vo = [ -141212  -989979 ];
%     %     case 'pasb2'
%     %         Vr = 277476; Vs =  [ -8.7237 3.38421 ]; Vo = [ 932306 -361672 ];
%     %     case 'ze11u' %='pasb2'
%     %         Vr = 277476; Vs =  [ -8.7237 3.38421 ]; Vo = [ 932306 -361672 ];
%     %     case 'ycWC0'
%     %         Vr = 291720; Vs =  [ 3.00029 -4.8053 ]; Vo = [ -529615 848238 ];
%     %     case 'lsvMi' %='ycWC0'
%     %         Vr = 291720; Vs =  [ 3.00029 -4.8053 ]; Vo = [ -529615 848238 ];
%     %     case '13215' %='ycWC0'
%     %         Vr = 291720; Vs =  [ 3.00029 -4.8053 ]; Vo = [ -529615 848238 ];
%     %     otherwise
%     %         warning('initial configuration unknown')
%     % end
%     
%     % Vt(1,:) = Vo(1)+Vs(1)*[1:10:nb_time*10];
%     % Vt(2,:) = Vo(2)+Vs(2)*[1:10:nb_time*10];
% 
%     % size_IP = 1e5;
%     % dist_to_center = sqrt(Vo*Vo');
%     
%     % inside1 = floor((dist_to_center-size_IP/2)/sqrt(Vs*Vs')/10);
%     % inside2 = floor(dist_to_center/sqrt(Vs*Vs')/10); 
%     % inside3 = floor((dist_to_center+size_IP/2)/sqrt(Vs*Vs')/10);
%     
%     % go_enter1 = floor(abs(dist_to_center-Vr-size_IP/2)/sqrt(Vs*Vs')/10);
%     % go_enter2 = floor(abs(dist_to_center-Vr)/sqrt(Vs*Vs')/10);
%     % go_enter3 = floor(abs(dist_to_center-Vr+size_IP/2)/sqrt(Vs*Vs')/10);
%     
%     % go_out1 = floor(abs(dist_to_center+Vr-size_IP/2)/sqrt(Vs*Vs')/10);
%     % go_out2 = floor(abs(dist_to_center+Vr)/sqrt(Vs*Vs')/10); 
%     % go_out3 = floor(abs(dist_to_center+Vr+size_IP/2)/sqrt(Vs*Vs')/10);
%     
%     % time_idx = [inside1;inside2;inside3;go_enter1;go_enter2;go_enter3;...
%     %     go_out1;go_out2;go_out3];
%     
%     % 1/ regimes: ballistic and browian
%     kfigs = 1;
%     figure(kfigs)
% 
%     loglog(10*delta_t,Diffu,'b');
%     xlabel('time (s)','FontSize',12);
%     ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex','FontSize',12);
%     hold on
%     grid on
% 
%     kfigs = kfigs+1;
%     figure(kfigs)
%     loglog(10*delta_t,Diffu,'b');
%     hold on
%     grid on
%     ax = gca;
%     ax.XLim =[1e4 5e5];
%     ax.YLim =[1e3 1e6];
%     ax.XTick = [12*3600:12*3600:3*24*3600 4*24*3600 5*24*3600];
%     ax.XTickLabel = {'1/2','1','1.5','2','2.5','3','4','5'};
%     xlabel('time (day)','FontSize',12);
%     ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex','FontSize',12);
% 
%     % time1 = 12*360; %(12h)
%     % loglog(time1:nb_time,Vd(time1:nb_time),'b');%'linewidth',2,'r');
%         
%     % hold on
%     % grid on
%     
%     % loglog(time1:nb_time,Vd_min(time1:nb_time),'--r');
%     % loglog(time1:nb_time,Vd_max(time1:nb_time),':r');
% 
%     % plot(inside1,Vd(inside2),'xb')
%     % plot(inside2,Vd(inside2),'xr')
%     % plot(inside3,Vd(inside2),'xb')
%     
%     % plot(go_enter1,Vd(go_enter2),'xb')
%     % plot(go_enter2,Vd(go_enter2),'xr')
%     % plot(go_enter3,Vd(go_enter2),'xb')
%     
%     % plot(go_out1,Vd(go_out2),'xb')
%     % plot(go_out2,Vd(go_out2),'xr')
%     % plot(go_out3,Vd(go_out2),'xb')
%     
%     % plot(2*inside2,Vd(2*inside2),'xb');
%     % ax = gca;
%     % ax.XLim = [8*360 8*24*360];
%     % ax.XTick = [12*360 24*360:24*360:8*24*360];
%     % ax.XTickLabel = {'1/2','1','2','3','4','5','6','7','8'};
%     % xlabel('time (day)');
%     % ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex');
%     % ax.XLabel.FontSize = 12;
%     % ax.YLabel.FontSize = 12;
%     
%     % kfigs = kfigs+1;    
%     % figure(kfigs)
%     % ax = gca;
%     % plot(Vt(1,inside2-400:inside2+400),Vt(2,inside2-400:inside2+400),'r')
%     % hold on
%     % grid on
%     % plot(Traj_MC(1,:),Traj_MC(2,:),'b')
%     % plot(Traj_MC(1,inside2),Traj_MC(2,inside2),'xr');
%     % plot(Vt(1,inside2-100),Vt(2,inside2-100),'xg');
%     % ax.XLim = [-40000 40000];
%     % ax.YLim = ax.XLim;
%     % ax.XTick = [-40000:10000:40000];
%     % ax.XTickLabel = {'-40','-30','-20','-10','0','10','20','30','40'};
%     % xlabel('distance (km)')
%     % ax.YTick = [-40000:10000:40000];
%     % ax.YTickLabel = {'-40','-30','-20','-10','0','10','20','30','40'};
%     % ylabel('distance (km)')
%     % ax.XLabel.FontSize = 12;
%     % ax.YLabel.FontSize = 12;
%     
%     % % large scale vortex and pack displacement:
%     % kfigs = kfigs+1;
%     % figure(kfigs)
%     % ax = gca;
%     % theta = 0:pi/32:2*pi;
%     % box = 10000*[-5 -5;5 -5;5 5;-5 5;-5 -5];
%     
%     % plot(Traj_MC(1,:),Traj_MC(2,:),'b')
%     % grid on
%     % hold on
%     
%     % plot(Vt(1,inside2)+Vr*cos(theta),Vt(2,inside2)+Vr*sin(theta),'g')
%     % plot(Vt(1,inside1),Vt(2,inside1),'xc')
%     % plot(Vt(1,inside2),Vt(2,inside2),'xc')
%     % plot(Vt(1,inside3),Vt(2,inside3),'xc')
%     
%     % plot(Traj_MC(1,inside2),Traj_MC(2,inside2),'xg')
%     % plot(Traj_MC(1,go_out2),Traj_MC(2,go_out2),'xk')
%     % plot(Traj_MC(1,go_enter2),Traj_MC(2,go_enter2),'xr')
%     % plot(Traj_MC(1,inside2)+box(:,1),Traj_MC(2,inside2)+box(:,2),'g')
%     % plot(Traj_MC(1,go_out2)+box(:,1),Traj_MC(2,go_out2)+box(:,2),'k')
%     % plot(Traj_MC(1,go_enter2)+box(:,1),Traj_MC(2,go_enter2)+box(:,2),'r')
%     
%     % plot(Vt(1,go_enter1)+Vr*cos(theta),Vt(2,go_enter1)+Vr*sin(theta),'b')
%     % plot(Vt(1,go_enter2)+Vr*cos(theta),Vt(2,go_enter2)+Vr*sin(theta),'r')
%     % plot(Vt(1,go_enter3)+Vr*cos(theta),Vt(2,go_enter3)+Vr*sin(theta),'b')
%     
%     % plot(Vt(1,go_out1)+Vr*cos(theta),Vt(2,go_out1)+Vr*sin(theta),'m')
%     % plot(Vt(1,go_out2)+Vr*cos(theta),Vt(2,go_out2)+Vr*sin(theta),'k')
%     % plot(Vt(1,go_out3)+Vr*cos(theta),Vt(2,go_out3)+Vr*sin(theta),'m')
%     
%     % ax.XLim = [-Vr-10000 Vr+10000];
%     % ax.YLim = [-Vr-10000 Vr+10000];
%     
%     % 2/ autocorrelation:
%     kfigs = kfigs+1;
%     figure(kfigs)
%     ax = gca;
%     plot(delta_t_a,ACF,'b');
%     grid on
%     hold on
%     ax.XLim = [0 5.5*24*360];
%     ax.XTick = [0:12*360:5.5*24*360];
%     ax.XTickLabel = {'0','1/2','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5'};
%     xlabel('time lag $\tau$ (day)','FontSize',12,'Interpreter','latex');
%     ylabel('normalized autocorrelation function','FontSize',12);
%     ax.XLim = [0 (simu_left_a/10+12*360)];
%     
%     leg = legend;
%     leg.String = {'$\tau \in \left[0,5.5j \right]$','$\tau \in \left[0,4.5j \right]$','$\tau \in \left[0,3j \right]$','$\tau \in \left[0,1.5j \right]$'};
%     leg.String = {'$60\%, a=2.5$','$70\%, a=1.5$','$70\%, a=2.5$','$80\%, a=1.5$'};
%     leg.Interpreter = 'latex';
%     title('Autocorrelation for 1 simu','FontSize',14);
%     leg.FontSize = 12;
%     
%     % 2bis/ Power Spectrum Density:
%     kfigs = kfigs+1;
%     figure(kfigs)
%     %N=6;
%     loglog(f,PSD,'b')
%     xlabel(' frequence (cyles per day) ','FontSize',12);
%     ylabel(' Power Spectral Density ','FontSize',12);
%     title('avg on 5 simus, ...','FontSize',12)
%     grid on;
%     hold on;
%     
%     % 3/ pair trajectory dispersion:
% %     kfigs = kfigs+1;
% %     figure(kfigs)
% %     ax = gca;
% %     loglog(time1:nb_time-1,disp_r(time1:nb_time-1),'b');
% %         
% %     hold on
% %     grid on
% %     
% %     plot(inside1,disp_r(inside1),'xb')
% %     plot(inside2,disp_r(inside2),'xr')
% %     plot(inside3,disp_r(inside3),'xb')
% %     
% %     plot(go_enter1,disp_r(go_enter1),'xb')
% %     plot(go_enter2,disp_r(go_enter2),'xr')
% %     plot(go_enter3,disp_r(go_enter3),'xb')
% %     
% %     plot(go_out1,disp_r(go_out1),'xb')
% %     plot(go_out2,disp_r(go_out2),'xr')
% %     plot(go_out3,disp_r(go_out3),'xb')
% %     
% %     plot(2*inside2,disp_r(2*inside2),'xb');
% %     
% %     loglog(time1:nb_time-1,disp_r_min_d(time1:nb_time-1),'--b');
% %     loglog(time1:nb_time-1,disp_r_max_d(time1:nb_time-1),':b');
% %     loglog(time1:nb_time-1,disp_r_min_S(time1:nb_time-1),'b--o','markersize',0.2);
% %     loglog(time1:nb_time-1,disp_r_max_S(time1:nb_time-1),'-.b');
% %     
% %     ax.XLim = [8*360 8*24*360];
% %     ax.XTick = [12*360 24*360:24*360:8*24*360];
% %     ax.XTickLabel = {'1/2','1','2','3','4','5','6','7','8'};
% %     xlabel('time (day)');
% %     ylabel('Pair trajectory dispersion');
% %     ax.XLabel.FontSize = 12;
% %     ax.YLabel.FontSize = 12;
%     
%     % 4/ distribution of tangential and normal acceleration
% %     kfigs = kfigs+1;
% %     figure(kfigs)
% %     ax = gca;
% %     tau_dis = [10 50 100 250:250:1000 2000:1000:5000 10000:5000:30000 40000 50000 75000 100000]/10;
% %     std_t = zeros(size(tau_dis));
% %     for i=1:size(tang_a,1)
% %         std_t(i) = std(tang_a(i,1:nb_time-tau_dis(i)-1));
% %         tang_std = tang_a(i,1:nb_time-tau_dis(i)-1)./std_t(i);
% %         xx = linspace(-7,7,500);
% %         [f_tang,x_tang] = ksdensity(tang_std,xx);
% % %         [f_norm,x_norm] = ksdensity(norm_a(i,1:nb_time-tau(i)-1));
% %         semilogy(x_tang,f_tang)
%         
% %         waitforbuttonpress
% %         if i==1
% %             hold on
% %             grid on
% %             ax.YLim = [1e-3 30];
% %             ax.XLim = [-4 4];
% %         end
% % %         plot(x_norm,f_norm,'--')
% %     end
%     
%     kfigs = kfigs+1;
%     figure(kfigs)
%     loglog(10*delta_t,acc(:,1),'b')
%     hold on
%     grid on
%     for j=2:length(q)
%         loglog(10*delta_t,acc(:,j),'b')
%     end
%     ylabel('$\left< \left( \frac{\Delta u}{\tau} \right)^q \right>$','Interpreter','latex','FontSize',12);
%     xlabel('time (s)','FontSize',12);
% %     loglog(10*tau_dis,std_t,'x-b');
% %     ax = gca;
% % %     ax.XTick = tau_dis;
% % %     ax.XTickLabel = {10 50 100 250:250:1000 2000:1000:5000 10000:5000:30000 40000 50000 75000 100000};
% %     ax.XLabel.String = 'time lag (s)';
% %     ax.YLabel.String = 'standard deviation';
% %     ax.XLabel.FontSize = 12;
% %     ax.YLabel.FontSize = 12;
% %     
% %     for i=1:size(tang_a,1)
% %         kfigs = kfigs+1;
% %         figure(kfigs)
% %         
% %         tang_a_p = tang_a(tang_a(i,1:nb_time-tau_dis(i)-1)>=0);
% %         tang_a_m = tang_a(tang_a(i,1:nb_time-tau_dis(i)-1)<=0);
% %         histogram(tang_a_p,'EdgeColor','r')
% %         ax = gca;
% %         ax.XScale = 'log';
% %         ax.YLim = [0 2500];
% %         hold on
% %         histogram(abs(tang_a_m),'EdgeColor','b')
% %         texttitle = strcat('tangential acceleration (',int2str(i),')');
% %         ax.Title.String = texttitle;
% %     end
%     
    
    %% saving:
    fprintf('saving phase in %s\n',mat_save_name)
    try
        save(mat_save_name,'nb_floes','nb_time','floe_traj_x','floe_traj_y',...
            'T','Mass_C','Area','Diffu','Diffu_min_s','Diffu_max_s','KE',...
            'Disper','ACF','ACF_min_s','ACF_max_s','PSD','f','acc',...
            'temporal_window','temporal_window_a')
    catch
        save(mat_save_name,'nb_floes','nb_time','floe_traj_x','floe_traj_y',...
            'T','Mass_C','Area','Diffu','Diffu_min_s','Diffu_max_s',...
            'Disper','ACF','ACF_min_s','ACF_max_s','PSD','f','acc',...
            'temporal_window','temporal_window_a')
    end

    % save(mat_save_name,'nb_floes','nb_time','floe_traj_x','floe_traj_y',...
    %     'idx_min_s','idx_max_s','Vd','Vd_min','Vd_max','time_idx',...
    %     'Vt','Traj_MC','tau_auto','KE','T','Mass_C','Area','Diffu',...
    %     'tau_dis','VS','ACF','f','Power','disp_r','idx_min_d','idx_max_d',...
    %     'disp_r_min_d','disp_r_max_d','disp_r_min_S','disp_r_max_S',...
    %     'idx_min_pS','idx_max_pS','tang_a','norm_a','std_t')
    
    % coucou = 10;
    % save(mat_save_name,'Vd_min','Vd_max','-append')
end 