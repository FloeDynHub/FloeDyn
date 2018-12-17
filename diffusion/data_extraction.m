%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data extraction (M. Rabatel IR 11-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of filenames must be a vector of string such as: 'simu_id'
%
% reminder: output filename are in the form: 'out_partial_simu_id.h5'
function data_extraction

% filenames = {'9VTBp','Ms3Dl','bSNvI','pasb2','ycWC0',...    % 60%, a=1.5
%     'T2PX8','c4AXm','n119i','7Tjnc','kcBWK',...  % 60, a=2.5
%     '3sbVt','UkdVH','ik9dr','oq73K','lo4si',...             % 70%, a=1.5
%     'Soq70','JlI88','iEmZN','Ihv7O','Q1C6I',...             % 70%, a=2.5
%     'gK8MP','Vg4W8','TTxQT','Ah4Rx','ITCTm'};               % 80%, a=1.5

% filenames = {'3sbVt','UkdVH','ik9dr','oq73K','lo4si','Soq70','JlI88','iEmZN','Ihv7O','Q1C6I'};

% filenames = {'gK8MP','Vg4W8','ITCTm'};

% data for statistics:
% temporal_window = 259200; % 3 days
% simu_left = 604800-temporal_window;
% delta_t = [10:10:1e2 150:50:1e3 1500:500:1e4 12500:5000:simu_left]/10;
% Col = {'b','r','k','g','m','c','y','--b','--r'};
% assert(length(Col) == length(filenames));

filenames = {'500f_15e-1a_9VTBp','500f_15e-1a_Ms3Dl',... % 500f, 10%, a=1.5
    '500f_15e-1a_bSNvI','500f_15e-1a_pasb2','500f_15e-1a_ycWC0'};%,...
%     '500f_25e-1a_9VTBp','500f_25e-1a_Ms3Dl',...   % 500f, 10%, a=2.5
%     '500f_25e-1a_bSNvI','500f_25e-1a_pasb2','500f_25e-1a_ycWC0'};            

wkey = 'null';

for idx_f=1:length(filenames)
    fprintf('Extraction of %s...\n',filenames{idx_f});
%     fileN = strcat('out_partial_',filenames{idx_f},'.h5');
    fileN = strcat('traj_analyse_',filenames{idx_f},'.mat');
%     traj_analyze(fileN);
    load(fileN);     wkey = 'gathering';

    if strcmp(wkey,'gathering')

        if idx_f==1
            meanDiffu = Diffu;
            meanDiffu_min_s = Diffu_min_s; meanDiffu_max_s = Diffu_max_s;
            
            meanACF = ACF; meanPSD = PSD;
            meanACF_min_s = ACF_min_s; meanACF_max_s = ACF_max_s;
            
            meanacc = acc;
        else
            meanDiffu = meanDiffu + Diffu;
            meanDiffu_min_s = meanDiffu_min_s + Diffu_min_s;
            meanDiffu_max_s = meanDiffu_max_s + Diffu_max_s;
            
            meanACF = meanACF+ACF; meanPSD = meanPSD+PSD;
            meanACF_min_s = meanACF_min_s + ACF_min_s;
            meanACF_max_s = meanACF_max_s + ACF_max_s;
            
            meanacc = meanacc+acc;
        end
    end
end

if strcmp(wkey,'gathering')

    Diffu       = meanDiffu/length(savemat);
    Diffu_min_s = meanDiffu_min_s/length(savemat);
    Diffu_max_s = meanDiffu_max_s/length(savemat);

    ACF         = meanACF/length(savemat);
    ACF_min_s   = meanACF_min_s/length(savemat);
    ACF_max_s   = meanACF_max_s/length(savemat);
    
    PSD         = meanPSD/length(savemat);

    ACC         = meanacc/length(savemat);
end

end