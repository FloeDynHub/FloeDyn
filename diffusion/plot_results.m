%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting results (M. Rabatel IR 11-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of file_name must be: 'file_name.h5'
%
function plot_results


for i=1:length(savemat)
    loadfile = strcat('traj_analyse_',savemat{i},'.mat');
    fprintf('loading results from %s...\n',savemat{i});
    load(loadfile);
    if i==1
        meanDiffu = Diffu;
        meanDiffu_min_s = Diffu_min_s; meanDiffu_max_s = Diffu_max_s;
        
        meanACF = ACF; 
        meanACF_min_s = ACF_min_s; meanACF_max_s = ACF_max_s;
        
        meanPSD = PSD;
        
        meanacc = acc; 
    else
        meanDiffu = meanDiffu + Diffu;
        meanDiffu_min_s = meanDiffu_min_s + Diffu_min_s;
        meanDiffu_max_s = meanDiffu_max_s + Diffu_max_s;
        
        meanACF = meanACF+ACF; 
        meanACF_min_s = meanACF_min_s + ACF_min_s;
        meanACF_max_s = meanACF_max_s + ACF_max_s;
        
        meanPSD = meanPSD+PSD;
        
        meanacc = meanacc+acc;
    end
end

Diffu       = meanDiffu/length(savemat);
Diffu_min_s = meanDiffu_min_s/length(savemat);
Diffu_max_s = meanDiffu_max_s/length(savemat);

ACF         = meanACF/length(savemat);
ACF_min_s   = meanACF_min_s/length(savemat);
ACF_max_s   = meanACF_max_s/length(savemat);
PSD         = meanPSD/length(savemat);

ACC         = meanacc/length(savemat);

save('traj_analyse_60p_25e-1a','temporal_window','temporal_window_a','ACF','ACF_min_s','ACF_max_s','Diffu','Diffu_min_s','Diffu_max_s','PSD')

% save('traj_analyse_60p_25e-1a.mat','ACF','PSD','ACC','-append');



time1 = 12*360;   %(12h)
N=6;

kfigs = 1;
Col = {'b','r','k','g','m'};
Colx = {'x-b','x-r','x-k','x-g','x-m'};
savemat = {'9VTBp','Ms3Dl','bSNvI','pasb2','ycWC0'}; % c=60, a=?
savemat = {'3sbVt','UkdVH','ik9dr','oq73K','lo4si'}; % c=70, a=1.5
savemat = {'Soq70','JlI88','iEmZN','Ihv7O','Q1C6I'}; % c=70, a=2.5
savemat = {'gK8MP','Vg4W8','TTxQT','Ah4Rx','ITCTm'}; % c=80, a=1.5

savemat = {'9VTBp','Ms3Dl','bSNvI','pasb2','ycWC0','3sbVt','UkdVH','ik9dr','ze11u'};

savemat = {'9VTBp','Ms3Dl','bSNvI','pasb2','ycWC0','3sbVt','UkdVH','ik9dr','oq73K','lo4si','Soq70','JlI88','iEmZN','Ihv7O','Q1C6I','TTxQT','Ah4Rx'};

%% individual plot
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:length(savemat)
        loadfile = strcat('traj_analyse_',savemat{i},'.mat');
        fprintf('loading results from %s...\n',savemat{i});
        load(loadfile);
        if i==1
            meanDiffu = Diffu; 
            meanDiffu_min_s = Diffu_min_s; meanDiffu_max_s = Diffu_max_s;
            
            meanACF = ACF; meanPSD = PSD;
            meanACF_min_s = ACF_min_s; meanACF_max_s = ACF_max_s;
            
            m_disp_r = disp_r; 
            m_disp_r_min_d = disp_r_min_d; 
            m_disp_r_max_d = disp_r_max_d;
            m_disp_r_min_S = disp_r_min_S;
            m_disp_r_max_S = disp_r_max_S;
          
            m_std_t = std_t;
            
            m_tang_a = tang_a;
            m_norm_a = norm_a;
            
            tau_dis = [10 50 100 250:250:1000 2000:1000:5000 10000:5000:30000 40000 50000 75000 100000]/10;
            std_t = std(tang_a(i,1:nb_time-tau_dis(i)-1));
            tang_std = tang_a(i,1:nb_time-tau_dis(i)-1)./std_t;
            xx = linspace(-7,7,500);
            [f_tang,~] = ksdensity(tang_std,xx);
            m_f_tang = f_tang;
        else
            meanDiffu = meanDiffu + Diffu; 
            meanDiffu_min_s = meanDiffu_min_s + Diffu_min_s; 
            meanDiffu_max_s = meanDiffu_max_s + Diffu_max_s;
            
            meanACF = meanACF+ACF; meanPSD = meanPSD+PSD;
            meanACF_min_s = meanACF_min_s + ACF_min_s; 
            meanACF_max_s = meanACF_max_s + ACF_max_s;
            
            m_disp_r = m_disp_r+disp_r; 
            m_disp_r_min_d = m_disp_r_min_d+disp_r_min_d;
            m_disp_r_max_d = m_disp_r_max_d+disp_r_max_d;
            m_disp_r_min_S = m_disp_r_min_S+disp_r_min_S;
            m_disp_r_max_S = m_disp_r_max_S+disp_r_max_S;
            
            m_std_t = m_std_t+std_t;
            
            m_tang_a = m_tang_a+tang_a;
            m_norm_a = m_norm_a+norm_a;
            
            tau_dis = [10 50 100 250:250:1000 2000:1000:5000 10000:5000:30000 40000 50000 75000 100000]/10;
            std_t = std(tang_a(i,1:nb_time-tau_dis(i)-1));
            tang_std = tang_a(i,1:nb_time-tau_dis(i)-1)./std_t;
            xx = linspace(-7,7,500);
            [f_tang,~] = ksdensity(tang_std,xx);
            m_f_tang = m_f_tang+f_tang;
        end

        %%% figure 1/ regimes: ballistic and browian
        figure(kfigs)
        loglog(time1:nb_time,Vd(time1:nb_time),Col{i});
        if i==1, hold on; grid on; end
        
%         loglog(time1:nb_time,Vd_min(time1:nb_time),'--r');
%         loglog(time1:nb_time,Vd_max(time1:nb_time),':r');
        
        plot(time_idx([1 3 4 6 7 9]),Vd(time_idx([1 3 4 6 7 9])),'xb');
        plot(time_idx([2 5 8]),Vd(time_idx([2 5 8])),'xr');
        plot(2*time_idx(5),Vd(2*time_idx(5)),'xb');
        
        if i==1
            ax = gca;
            ax.XLim = [8*360 8*24*360];
            ax.XTick = [12*360 24*360:24*360:8*24*360];
            ax.XTickLabel = {'1/2','1','2','3','4','5','6','7','8'};
            xlabel('time (day)');
            ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex');
            ax.XLabel.FontSize = 12;
            ax.YLabel.FontSize = 12;
        end
        
        %%% figure % 2/ autocorrelation:
        kfigs = kfigs+1;
        figure(kfigs)
        plot(tau_auto,ACF,Col{i});
        if i==1, hold on; grid on; end
        if i==1
            ax = gca;
            ax.XLim = [0 5.5*24*360];
            ax.XTick = [0:12*360:5.5*24*360];
            ax.XTickLabel = {'0','1/2','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5'};
            xlabel('time lag (day)');
            ylabel('normalized autocorrelation function');
            ax.XLabel.FontSize = 12;
            ax.YLabel.FontSize = 12;
        end
    
        %%% figure % 2bis/ Power Spectrum Density:
        kfigs = kfigs+1;
        figure(kfigs)
        loglog(f,Power(1:2^(N-1)),Colx{i})
        if i==1, hold on; grid on; end
        if i==1
            xlabel(' Frequency (Hz) ')
            ylabel(' Magnitude (w) '),
            title('Power Spectral Density')
        end
        
        %%% figure % 3/ pair trajectory dispersion:
        kfigs = kfigs+1;
        figure(kfigs)
        loglog(time1:nb_time,disp_r(time1:nb_time),Col{i});
        if i==1, hold on; grid on; end
        
        plot(time_idx([1 3 4 6 7 9]),Vd(time_idx([1 3 4 6 7 9])),'xb');
        plot(time_idx([2 5 8]),Vd(time_idx([2 5 8])),'xr');
        plot(2*time_idx(5),Vd(2*time_idx(5)),'xb');
        
%         loglog(time1:nb_time,disp_r_min_d(time1:nb_time),'--b');
%         loglog(time1:nb_time,disp_r_max_d(time1:nb_time),':b');
%         loglog(time1:nb_time,disp_r_min_S(time1:nb_time),'b--o','markersize',0.2);
%         loglog(time1:nb_time,disp_r_max_S(time1:nb_time),'-.b');
        
        if i==1
            ax = gca;
            ax.XLim = [8*360 8*24*360];
            ax.XTick = [12*360 24*360:24*360:8*24*360];
            ax.XTickLabel = {'1/2','1','2','3','4','5','6','7','8'};
            xlabel('time (day)');
            ylabel('Pair trajectory dispersion');
            ax.XLabel.FontSize = 12;
            ax.YLabel.FontSize = 12;
        end

        %%% figure % 4/ distribution of tangential and normal acceleration
        kfigs = kfigs+1;
        figure(kfigs)
        loglog(10*tau_dis,std_t,Colx{i});
        if i==1, hold on; grid on; end
        if i==1
            ax = gca;
            ax.XLabel.String = 'time lag (s)';
            ax.YLabel.String = 'standard deviation';
            ax.XLabel.FontSize = 12;
            ax.YLabel.FontSize = 12;
        end
    
        kfigs = 1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('traj_analyse_9VTBp.mat')
%      Vr = 374848; Vs =  [ -9.10191 0.813075 ]; Vo = [ 996034  -88975.9 ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('traj_analyse_Ms3Dl.mat')
%     Vr = 325639; Vs =  [ 1.26598 8.87528 ]; Vo = [ -141212  -989979 ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('traj_analyse_bSNvI.mat')
%     Vr = 342577; Vs =  [ 3.73322 9.1694 ]; Vo = [ -377083 -926180 ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('traj_analyse_pasb2.mat')
%     Vr = 277476; Vs =  [ -8.7237 3.38421 ]; Vo = [ 932306 -361672 ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     load('traj_analyse_ycWC0.mat')
%     Vr = 291720; Vs =  [ 3.00029 -4.8053 ]; Vo = [ -529615 848238 ];

%% multiple plots:

kfigs = kfigs + 20;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Diffu       = meanDiffu/length(savemat);
    Diffu_min_s = meanDiffu_min_s/length(savemat);
    Diffu_max_s = meanDiffu_max_s/length(savemat);

    ACF         = meanACF/length(savemat);
    ACF_min_s   = meanACF_min_s/length(savemat);
    ACF_max_s   = meanACF_max_s/length(savemat);
    PSD         = meanPSD/length(savemat);
    
    m_disp_r = m_disp_r/length(savemat);
    m_disp_r_min_d = m_disp_r_min_d/length(savemat);
    m_disp_r_max_d = m_disp_r_max_d/length(savemat);
    m_disp_r_min_S = m_disp_r_min_S/length(savemat);
    m_disp_r_max_S = m_disp_r_max_S/length(savemat);
    
    m_std_t = m_std_t/length(savemat);
    
    m_tang_a = m_tang_a/length(savemat);
    m_norm_a = m_norm_a/length(savemat);
    
    m_f_tang = m_f_tang/length(savemat);
    
    %%% figure % 1/ regimes: ballistic and browian
    figure(kfigs)
    loglog(time1:nb_time,meanVd(time1:nb_time),'b');
    
    hold on
    grid on
    
    loglog(time1:nb_time,meanVd_min(time1:nb_time),'--r');
    loglog(time1:nb_time,meanVd_max(time1:nb_time),':r');

    ax = gca;
    ax.XLim = [8*360 8*24*360];
    ax.XTick = [12*360 24*360:24*360:8*24*360];
    ax.XTickLabel = {'1/2','1','2','3','4','5','6','7','8'};
    xlabel('time (day)');
    ylabel('$\left< {r^\prime}^2 \right>$','Interpreter','latex');
    title('average on 5 simus')
    ax.XLabel.FontSize = 12;
    ax.YLabel.FontSize = 12;
    
    %%% figure % 2/ autocorrelation:
    kfigs = kfigs+1;
    figure(kfigs)
    plot(tau_auto,meanACF,'b');
    grid on;
    
    ax = gca;
    ax.XLim = [0 5.5*24*360];
    ax.XTick = [0:12*360:5.5*24*360];
    ax.XTickLabel = {'0','1/2','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5'};
    xlabel('time lag (day)');
    ylabel('normalized autocorrelation function');
    title('average on 5 simus')
    ax.XLabel.FontSize = 12;
    ax.YLabel.FontSize = 12;
    
     %%% figure % 2bis/ Power Spectrum Density:
     kfigs = kfigs+1;
     figure(kfigs)
     %N=6;
     loglog(f,meanPower(1:2^(N-1)),'x-b')
     xlabel('  Frequency (Hz)')
     ylabel(' Magnitude (w)'),
     title('Power Spectral Density (average on 5 simus)')
     grid on;
    
     
     %%% figure % 3/ pair trajectory dispersion:
     kfigs = kfigs+1;
     figure(kfigs)
     loglog(time1:nb_time,m_disp_r(time1:nb_time),'b');
     
     hold on
     grid on
     
     loglog(time1:nb_time,m_disp_r_min_d(time1:nb_time),'--b');
     loglog(time1:nb_time,m_disp_r_max_d(time1:nb_time),':b');
     loglog(time1:nb_time,m_disp_r_min_S(time1:nb_time),'b--o','markersize',0.2);
     loglog(time1:nb_time,m_disp_r_max_S(time1:nb_time),'-.b');
     
     ax.XLim = [8*360 8*24*360];
     ax.XTick = [12*360 24*360:24*360:8*24*360];
     ax.XTickLabel = {'1/2','1','2','3','4','5','6','7','8'};
     xlabel('time (day)');
     ylabel('Pair trajectory dispersion');
     title('average on 5 simus');
     ax.XLabel.FontSize = 12;
     ax.YLabel.FontSize = 12;
     
     
     %%% figure % 4/ distribution of tangential and normal acceleration
     kfigs = kfigs+1;
     figure(kfigs)
     loglog(10*tau_dis,m_std_t,'x-b');
     ax = gca;
     ax.XLabel.String = 'time lag (s)';
     ax.YLabel.String = 'standard deviation';
     title('average on 5 simus');
     ax.XLabel.FontSize = 12;
     ax.YLabel.FontSize = 12;
     
     
     kfigs = kfigs+1;
     figure(kfigs)
     ax = gca;
     tau_dis2 = [250 1000 5000 50000 100000]/10;
     cc = [4 7 11 18 20];
     for i=1:length(cc)
         std_t = std(m_tang_a(cc(i),1:nb_time-tau_dis2(i)-1));
         tang_std = m_tang_a(cc(i),1:nb_time-tau_dis2(i)-1)./std_t;
         xx = linspace(-7,7,500);
         [f_tang,x_tang] = ksdensity(tang_std,xx);
         %         [f_norm,x_norm] = ksdensity(norm_a(i,1:nb_time-tau(i)-1));
         semilogy(x_tang,f_tang)
         
         waitforbuttonpress
         if i==1
             hold on
             grid on
             ax.YLim = [1e-3 30];
             ax.XLim = [-4 4];
         end
         %         plot(x_norm,f_norm,'--')
     end
     xlabel('$a_{\parallel} / \left< {a_{\parallel}}^2 \right>^{1/2}$','Interpreter','latex');
     leg = legend;
     leg.String = {'$\tau = 250s$','$\tau = 1000s$','$\tau = 5000s$','$\tau = 50000s$','$\tau = 100000s$'};
    leg.Interpreter = 'latex';
    leg.FontSize = 12;
    title('PDF of tangential accelerations');
end