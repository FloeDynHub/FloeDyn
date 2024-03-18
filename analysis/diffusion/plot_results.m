%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting results (M. Rabatel IR 11-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of file_name must be: 'file_name.h5'
%
function plot_results

% filenames = {'60p_a_unknown','60p_15e-1a','60p_25e-1a','70p_15e-1a','70p_25e-1a',...
%     '80p_15e-1a'};
filenames = {'60p_a_unknown','60p_25e-1a','70p_15e-1a','70p_25e-1a',...
    '80p_15e-1a'};

Col     = {'m','b','r','k','g'};
Leg_txt = {'60%, $\alpha=1.5$','60%, $\alpha=2.5$','70%, $\alpha=1.5$','70%, $\alpha=2.5$',...
    '80%, $\alpha=1.5$'};

for idx_f=1:length(filenames)
    loadfile = strcat('traj_analyse_',filenames{idx_f},'.mat');
    fprintf('loading results from %s...\n',filenames{idx_f});
    load(loadfile);
    
    temporal_window = 259200; % 3 days
    simu_left = 604800-temporal_window;
    % origin time along the trajectory:
    t_o = floor(linspace(1,(temporal_window/10),500)); % divided by 10 (min time_step)
    % time step:
    delta_t = [10:10:1e2 150:50:1e3 1500:500:1e4 12500:5000:simu_left]/10;
    
    temporal_window_a = 345600; % 4 days
    simu_left_a = 604800-temporal_window_a;
    % origin time along the trajectory:
    t_o_a = floor(linspace(1,(temporal_window_a/10),500)); % divided by 10 (min time_step)
    % time step:
    Frq = 48; % 48 per day -> 30 min -> 1800s
    day_in_s = 86400;
    delta_t_a = [0:day_in_s/Frq:simu_left_a]/10;
    
    %% Plot Phase
    % 1/ regimes: ballistic and browian
    fprintf('Plot phase ...\n');
    kfigs = 1;
    figure(kfigs)
    
    loglog(10*delta_t,Diffu,Col{idx_f});
    
    if idx_f == 1
        xlabel('time (s)','FontSize',12);
        ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex','FontSize',12);
        hold on
        grid on
    end
        
    kfigs = kfigs+1;
    figure(kfigs)
    
    loglog(10*delta_t,Diffu,Col{idx_f});
    
    if idx_f == 1
        hold on
        grid on
        ax = gca;
        ax.XLim =[1e4 5e5];
        ax.YLim =[1e3 1e6];
        ax.XTick = [12*3600:12*3600:3*24*3600 4*24*3600 5*24*3600];
        ax.XTickLabel = {'1/2','1','1.5','2','2.5','3','4','5'};
        xlabel('time (day)','FontSize',12);
        ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex','FontSize',12);
    end
    
    % 2/ autocorrelation:
    kfigs = kfigs+1;
    figure(kfigs)
    
    plot(delta_t_a,ACF,Col{idx_f});
    
    if idx_f == 1
        grid on
        hold on
        ax = gca;
        ax.XLim = [0 5.5*24*360];
        ax.XTick = [0:12*360:5.5*24*360];
        ax.XTickLabel = {'0','1/2','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5'};
        xlabel('time lag $\tau$ (day)','FontSize',12,'Interpreter','latex');
        ylabel('normalized autocorrelation function','FontSize',12);
        ax.XLim = [0 (simu_left_a/10+12*360)];
    end
    
    % 2bis/ Power Spectrum Density:
    kfigs = kfigs+1;
    figure(kfigs)

    L = length(delta_t_a);
    f = Frq*(1:ceil(L/2+1))/L;
    loglog(f,PSD,Col{idx_f})
    
    if idx_f == 1
        xlabel(' frequence (cyles per day) ','FontSize',12);
        ylabel(' Power Spectral Density ','FontSize',12);
        title('avg on 5 simus, ...','FontSize',12)
        grid on;
        hold on;
    end
    
    % 3/ acceleration
    kfigs = kfigs+1;
    figure(kfigs)
    
    q = [0.5 1 1.5 2 2.5 3];
    
    loglog(10*delta_t,ACC(:,1),Col{idx_f})
    
    if idx_f == 1
        hold on
        grid on
        ylabel('$\left< \left( \frac{\Delta u}{\tau} \right)^q \right>$','Interpreter','latex','FontSize',12);
        xlabel('time (s)','FontSize',12);
    end
    for j=2:length(q)
        loglog(10*delta_t,ACC(:,j),Col{idx_f})
    end
end

for i=1:kfigs-1
    figure(kfigs)
    leg = legend;
    leg.Interpreter = 'latex';
    leg.FontSize = 12;
    leg.String = Leg_txt;
end

end