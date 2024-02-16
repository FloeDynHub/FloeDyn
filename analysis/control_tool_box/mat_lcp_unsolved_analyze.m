%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Analyzing matrix from LCP unsolved (M. Rabatel IR 01-2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% the format of file_name must be: 'file_name.h5'
%
% Information on LCP(M,q)
%
% n_h5: numeration in h5 files
% n_mat: numeration in matlab files
%
% For matrix M from solved or unsolved LCP:
% n_mat = find(Idx_sol_alive==Ep(n_h5))
% n_mat = find(Idx_sol_alive_uns==Ep_uns(n_h5))
% 
% n_h5 = find(Ep==Idx_sol_alive(n_mat))
% n_h5 = find(Ep_uns==Idx_sol_alive_uns(n_mat))
% 
% Note: if find(Idx_sol_alive==Ep(n_h5)) is empty, then this matrix is
% similar to another one. To find it, do:
% KK = 0;
% for i=1:9812
%     if any(list_same_lcp{i}==Ep(n_h5))
%         KK=i;
%         break;
%     end
% end
%
function mat_lcp_unsolved_analyze(file_name)

    %% Preliminary infos:
    GG_1 = 'solved';
    GG_2 = 'unsolved';
    group_name_1 = 'M';
    group_name_2 = 'q';
    group_name_3 = 'z';
    Stats_solver = 'Which solver';
    Stats_err = 'LCP error';
    Loop_Info = 'Contact Graph Info';
        
    suf = split(file_name,'LCP_stats_');
    suffix = split(suf{2},'.h5');
    
    info = h5info(file_name);
    nb_LCP_solved = size(info.Groups(1).Groups(1).Datasets,1);
    nb_LCP_unsolved = size(info.Groups(2).Groups(1).Datasets,1);
    total_nb_lcp = nb_LCP_solved + nb_LCP_unsolved;
    
    fprintf('number of solved LCP: %d \n number of diff unsolved LCP: %d \n',...
        nb_LCP_solved, nb_LCP_unsolved);
    
    %% Storing in matlab format:
    fprintf('Storing phase...\n');
    M = cell(nb_LCP_solved,1); M_uns = cell(nb_LCP_unsolved,1);
    Q = M; Q_uns = M_uns; Z = M; Z_uns = M_uns;
    Ep = zeros(total_nb_lcp,1); Ep_uns = zeros(total_nb_lcp,1);

    j=1; k=1;
    for i=1:total_nb_lcp
        if mod(i,1000)==0; fprintf('%.4g%%\n',i/total_nb_lcp*100); end
        try
            member_name = strcat('/',GG_1,'/',group_name_1,'/',int2str(i));
            A = h5read(file_name,member_name); % read en transposed ??!!
            M{k} = A';
            member_name = strcat('/',GG_1,'/',group_name_2,'/',int2str(i));
            Q{k} = h5read(file_name,member_name); % read en transposed ??!! 
            member_name = strcat('/',GG_1,'/',group_name_3,'/',int2str(i));
            Z{k} = h5read(file_name,member_name); % read en transposed ??!!   
            Ep(i) = k;
            k = k+1;
        catch
        end
        
        try
            member_name = strcat('/',GG_2,'/',group_name_1,'/',int2str(i));
            A = h5read(file_name,member_name); % read en transposed ??!!
            M_uns{j} = A';
            member_name = strcat('/',GG_2,'/',group_name_2,'/',int2str(i));
            Q_uns{j} = h5read(file_name,member_name); % read en transposed ??!!
            member_name = strcat('/',GG_2,'/',group_name_3,'/',int2str(i));
            Z_uns{j} = h5read(file_name,member_name); % read en transposed ??!!
            Ep_uns(i) = j;
            j = j+1;
        catch
        end
    end
    
    try
        Contact_Graph = h5read(file_name,strcat('/',Loop_Info));
        Contact_Graph = Contact_Graph';
    catch
        Contact_Graph = [];
    end
    
    %% Computation of number of different LCP among the solved one:
    % solved one:
    ans_rem_sim_LCP = input('Do you want to (eventually) remove similar LCP? a string: [y/n]\n ');
    if strcmp(ans_rem_sim_LCP,'y')
    fprintf('Reduction of the solved LCP number by elimination of redundant LCP...\n');
    list_left_lpc = 1:nb_LCP_solved; list_same_lcp = cell(nb_LCP_solved,1);
    list_same_Q = [];
    
    while size(list_left_lpc,2)>0
        idx_lcp = list_left_lpc(1);
        if mod(idx_lcp,500)==0; fprintf('%.4g%%\n',idx_lcp/nb_LCP_solved*100); end
        
        lcp = M{idx_lcp};
        list_left_lpc = setdiff(list_left_lpc,idx_lcp);
        
        list_temp = [];
        for i=1:size(list_left_lpc,2)
            lcp_test = M{list_left_lpc(i)};
            
            if size(lcp,1) == size(lcp_test,1)
                DM1 = lcp(1:3*size(lcp,1)/4,1:3*size(lcp,1)/4);
                DM2 = lcp_test(1:3*size(lcp_test,1)/4,1:3*size(lcp_test,1)/4);
        
                Diff = DM1 - DM2;
                Diff_rel = Diff./min(abs(DM1),abs(DM2));
                Diff_rel_m = max(abs(Diff_rel),0);
                N_Diff_rel_m = norm(Diff_rel_m);
                
                if N_Diff_rel_m==0
                    list_temp = [list_temp;list_left_lpc(i)]; 
                    
                    vec_temp = Q{idx_lcp} - Q{list_left_lpc(i)};
                    N_sim_rel_v = norm(vec_temp);
                    if N_sim_rel_v==0
                        list_same_Q = [list_same_Q;idx_lcp list_left_lpc(i)];
                    end
                end
            end
        end
        list_same_lcp{idx_lcp} = list_temp;
        
        list_left_lpc = setdiff(list_left_lpc,list_temp);
    end
    
    fprintf('Be careful!! There are %d LCP with exactly same M and q!!\n\n',size(list_same_Q,2));
%     assert(size(list_same_Q,1)==0,'Exactly same LCP(q,M) is solved more than once!!');
    
    list_still_alive = 1:nb_LCP_solved; Idx_sol_alive = [];
    for i=1:nb_LCP_solved
        if ismember(i,list_still_alive)
            Idx_sol_alive = [Idx_sol_alive; i];
            ll = list_same_lcp{i};
            ll = [ll; i];
            list_still_alive = setdiff(list_still_alive,ll);
        end
    end
    
    M_s_alive = cell(size(Idx_sol_alive));
    Q_s_alive = cell(size(Idx_sol_alive));
    for i=1:size(Idx_sol_alive,1)
        M_s_alive{i} = M{Idx_sol_alive(i)};
        Q_s_alive{i} = Q{Idx_sol_alive(i)};
    end
    
    clear M Q
    M = M_s_alive; Q = Q_s_alive;
    clear M_s_alive Q_s_alive
    
    list_same_lcp_sol = list_same_lcp;
    fprintf('There are %d solved LCP with same M, yet there are %d with same M and q.\n'...
        ,nb_LCP_solved-size(Idx_sol_alive,1), size(list_same_Q,2));
    nb_LCP_solved = size(M,1);
    
    % unsolved one:
    fprintf('Reduction of the unsolved LCP number by elimination of redundant LCP...\n');
    list_left_lpc = 1:nb_LCP_unsolved; list_same_lcp = cell(nb_LCP_unsolved,1);
    
    while size(list_left_lpc,2)>0
        idx_lcp = list_left_lpc(1);
        if mod(idx_lcp,500)==0; fprintf('%.4g%%\n',idx_lcp/nb_LCP_unsolved*100); end
        
        lcp = M_uns{idx_lcp};
        list_left_lpc = setdiff(list_left_lpc,idx_lcp);
        
        list_temp = [];
        for i=1:size(list_left_lpc,2)
            lcp_test = M_uns{list_left_lpc(i)};
            
            if size(lcp,1) == size(lcp_test,1)
                DM1 = lcp(1:3*size(lcp,1)/4,1:3*size(lcp,1)/4);
                DM2 = lcp_test(1:3*size(lcp_test,1)/4,1:3*size(lcp_test,1)/4);
        
                Diff = DM1 - DM2;
                Diff_rel = Diff./min(abs(DM1),abs(DM2));
                Diff_rel_m = max(abs(Diff_rel),0);
                N_Diff_rel_m = norm(Diff_rel_m);
                
                if N_Diff_rel_m==0
                    vec_temp = Q_uns{idx_lcp} - Q_uns{list_left_lpc(i)};
                    N_sim_rel_v = norm(vec_temp);
                    if N_sim_rel_v==0
                        list_temp = [list_temp; list_left_lpc(i)];
                    end
                end
            end
        end
        list_same_lcp{idx_lcp} = list_temp;
        
        list_left_lpc = setdiff(list_left_lpc,list_temp);
    end
    
    list_still_alive = 1:nb_LCP_unsolved; Idx_sol_alive_uns = [];
    for i=1:nb_LCP_unsolved
        if ismember(i,list_still_alive)
            Idx_sol_alive_uns = [Idx_sol_alive_uns; i];
            ll = list_same_lcp{i};
            ll = [ll; i];
            list_still_alive = setdiff(list_still_alive,ll);
        end
    end
    
    M_s_alive = cell(size(Idx_sol_alive_uns));
    Q_s_alive = cell(size(Idx_sol_alive_uns));
    for i=1:size(Idx_sol_alive_uns,1)
        M_s_alive{i} = M_uns{Idx_sol_alive_uns(i)};
        Q_s_alive{i} = Q_uns{Idx_sol_alive_uns(i)};
    end
    
    clear M_uns Q_uns
    M_uns = M_s_alive; Q_uns = Q_s_alive;
    clear M_s_alive Q_s_alive
    
    list_same_lcp_unsol = list_same_lcp;
    fprintf('There are %d unsolved LCP with exactly same M and q.\n',...
        nb_LCP_unsolved-size(Idx_sol_alive_uns,1));
    nb_LCP_unsolved = size(M_uns,1);
    
    % recovery of solution (z) of LCP(q,M) coming from C++ implementation:
%     Z = cell(nb_LCP_solved,1); Z_uns = cell(nb_LCP_unsolved,1);
%     for i=1:nb_LCP_unsolved
%         n_h5 = find(Ep_uns==Idx_sol_alive_uns(i));
%         member_name = strcat('/',GG_2,'/',group_name_3,'/',int2str(n_h5));
%         Z_uns{i} = h5read(file_name,member_name); % read en transposed ??!!
%     end
%     
%     for i=1:nb_LCP_solved
%         n_h5 = find(Ep==Idx_sol_alive(i));
%         member_name = strcat('/',GG_1,'/',group_name_3,'/',int2str(n_h5));
%         Z{i} = h5read(file_name,member_name); % read en transposed ??!!
%     end
    
    %% Computation of number of different LCP:
    fprintf('Association between unsolved and solved LCP...\n');
    Idx_uns = []; k = 1;
    for i=1:nb_LCP_unsolved
        if mod(i,10)==0; fprintf('%.4g%%\n',i/nb_LCP_unsolved*100); end
        
        A = M_uns{i};
        % recovering the corresponding solved matrix:
        n_h5 = find(Ep_uns==Idx_sol_alive_uns(i));
        for j=n_h5+1:total_nb_lcp
            n_mat = find(Idx_sol_alive==Ep(j));
            if size(n_mat,1)==1
                B = M{n_mat};
        
                % Comparison between Delassus part:
                if size(A,1) == size(B,1)
                    DM_uns = A(1:3*size(A,1)/4,1:3*size(A,1)/4);
                    DM = B(1:3*size(B,1)/4,1:3*size(B,1)/4);

                    Diff = DM_uns - DM;
                    Diff_rel = Diff./min(abs(DM_uns),abs(DM));
                    Diff_rel_m = max(abs(Diff_rel),0);
                    N_Diff_rel_m = norm(Diff_rel_m);
                    if N_Diff_rel_m==0
                        vec_temp = Q_uns{i} - Q{n_mat};
                        N_sim_rel_v = norm(vec_temp);
                        Idx_uns(k,1) = N_sim_rel_v;
                        Idx_uns(k,2:3) = [i,n_mat];
                        k=k+1;
                    end
                end
            end
        end
    end
    
    %% Comparison with relative velocities q:
    fprintf('Separation with never solved LCP...\n');
    if size(Idx_uns,1)>0
        is_same_LCP = unique(Idx_uns(Idx_uns(:,1)==0,2));
    else
        is_same_LCP = [];
    end
    
    isnt_same_LCP = setdiff([1:nb_LCP_unsolved]',is_same_LCP);
    
    Idx_diff_uns_sol = zeros(size(isnt_same_LCP,1),4);
    if size(isnt_same_LCP,1) < 500
        for i=1:size(isnt_same_LCP,1)
            if mod(i,20)==0; fprintf('%.4g%%\n',i/nb_LCP_unsolved*100); end
            
            A = M_uns{isnt_same_LCP(i)};
            n_h5 = find(Ep_uns==Idx_sol_alive_uns(isnt_same_LCP(i)));
            min_N = inf;
            for j=n_h5+1:total_nb_lcp
                n_mat_s = find(Idx_sol_alive==Ep(j));
                if size(n_mat_s,1)==1
                    B = M{n_mat_s};
                    
                    % Comparison between Delassus part:
                    if size(A,1) == size(B,1)
                        DM_uns = A(1:3*size(A,1)/4,1:3*size(A,1)/4);
                        DM = B(1:3*size(B,1)/4,1:3*size(B,1)/4);
                        
                        Diff = DM_uns - DM;
                        Diff_rel = Diff./min(abs(DM_uns),abs(DM));
                        Diff_rel_m = max(abs(Diff_rel),0);
                        N_Diff_rel_m = norm(Diff_rel_m);
                        
                        if N_Diff_rel_m < min_N
                            min_N = N_Diff_rel_m;
                            Idx_M_sol_closest = n_mat_s;
                        end
                    end
                end
            end
            if min_N~=Inf
                vec_temp = Q_uns{isnt_same_LCP(i)} - Q{Idx_M_sol_closest};
                norm_Q = norm(vec_temp);
                k = find(Ep==Idx_sol_alive(Idx_M_sol_closest));
            else
                k = 0;
                norm_Q = -1;
            end
            norm_M = min_N;
            
            Idx_diff_uns_sol(i,:) = [norm_M norm_Q n_h5 k];
        end
    end
    end
    
    %% Recovering statistic on which solver and how many calls
    %|     1     |     2      |     3      |   4   |          5         |   6    |      8      |        9       |
    %|:---------:|:----------:|:----------:|:-----:|:------------------:|:------:|:-----------:|:--------------:|
    %| LCP error | nb attempt | nb perturb | nb SR | nb adj cone failed | lexico | idx failure | last technique |
    %
    % with: \e nb for number, \e SR for secondary ray, \e adj \e cone \e failed for the failure of the 
    % method consisting in going through an adjacent cone, \e lexico is true if the lexicographic 
    % ordering is used during at least one Lemke's algorithm, \e idx \e perturb for the index of the
    % matrix perturbation (see LCPSolver::matrix_perturbation(const std::size_t dim , matrix &M, const double alpha, const int Idx_perturb)),
    % \e idx \e failure for the source of the LCP error (see LCPSolver::which_failure( vector<real_type> Err, bool Is_pos_rel_norm_vel ))
    % and \e last \e technique for the last SR or perturbation used before solving or last attempt done.
    fprintf('Statistics on LCP solvers...\n');
    member_name = strcat('/',GG_2,'/',Stats_err);
    data_err_uns = h5read(file_name,member_name);
    member_name = strcat('/',GG_1,'/',Stats_err);
    data_err = h5read(file_name,member_name);
    
    data_err_uns = data_err_uns';    
    data_err = data_err';
    
    val_max_attempt = input('What is the maximum number of attempt for solving? \n');
    
    if nb_LCP_solved>0
        fprintf('the mean LCP error for solved LCP is: %.4e\n',mean(data_err(:,1)));
        for i=1:val_max_attempt
            fprintf('number of LCP solved after %d attempts: %d (%.4g%%)\n', i, ...
                sum(data_err(:,2)==i), sum(data_err(:,2)==i)/nb_LCP_solved*100);
        end
        
        fprintf('\n');
        for i=1:max(data_err(:,3))
            fprintf('number of LCP solved after %d perturbations: %d (%.4g%%) \n', ...
                i, size(find(data_err(:,3)==i),1), size(find(data_err(:,3)==i),1)/...
                nb_LCP_solved*100 );
        end
       
        fprintf('\n');
        for i=1:max(data_err(:,4))
            fprintf('number of LCP solved after %d Secondary Ray: %d (%.4g%%) \n', ...
                i, size(find(data_err(:,4)==i),1), size(find(data_err(:,4)==i),1)/...
                nb_LCP_solved*100 );
        end
        
        fprintf('\n');
        Using_perturb = find(data_err(:,3)>0);
        b_SR_PR = size(find(data_err(Using_perturb,4)>0),1);
        fprintf('number of LCP solved after using both SR and perturbations: %d (%.4g%%) \n', ...
            b_SR_PR, b_SR_PR/nb_LCP_solved*100 );
        
        fprintf('\n');
        fprintf('number of LCP solved although SR failed: %d (%.4g%%) \n',...
            size(find(data_err(:,5)>0),1), size(find(data_err(:,5)>0),1)/...
            nb_LCP_solved*100);
        
        fprintf('\n');
        fprintf('number of LCP solved with lexicographic ordering: %d (%.4g%%) \n',...
            size(find(data_err(:,6)>0),1), size(find(data_err(:,6)>0),1)/...
            nb_LCP_solved*100);
        
    end
    
    fprintf('\n');
    count_degenerate = 0; count_degenerate_FC = 0;
    for i=1:nb_LCP_solved
        W = M{i}*Z{i} + Q{i};
        supp_zer_sol = setdiff( 1:1:size(Z{i},1) , ...
            union( find(abs(Z{i})>1e-9), find(abs(W)>1e-9) ) );
        if size(supp_zer_sol,2)>0
            count_degenerate = count_degenerate+1;
        end
        
        % for FC formulation (without \alpha-basis and only one
        % \beta-basis)
        dim_FC = size(Z{i},1)/2; dim_nc = size(Z{i},1)/4;
        E_tmp = union( find(abs(Z{i})>1e-9),  find(abs(W)>1e-9) );
        E_tmp = E_tmp( E_tmp < 3*dim_nc+1 );
        E_tmp = E_tmp( E_tmp > dim_nc ); 
        rest_beta = setdiff( 1:1:dim_nc , ceil( (E_tmp-dim_nc) / 2 ) );
        rest_alpha = setdiff( 1:1:dim_nc , ...
            union( find( abs( Z{i}(1:dim_nc) ) > 1e-9 ), ...
            find( abs( W(1:dim_nc) ) > 1e-9 ) ) );
        
        supp_zer_sol = union( rest_alpha , rest_beta );
        if size(supp_zer_sol,2)>0
            count_degenerate_FC = count_degenerate_FC+1;
        end
    end
    
    if nb_LCP_solved==0
        percent_c_d = 100; percent_c_d_FC = 100;
    else
        percent_c_d     = count_degenerate      / nb_LCP_solved*100;
        percent_c_d_FC  = count_degenerate_FC   / nb_LCP_solved*100;
    end
    fprintf('%.4g%% (%d) of the solved LCP are degenerate (on the (APS) form) \n', ...
        percent_c_d, count_degenerate );
    fprintf('%.4g%% (%d) of the solved LCP are degenerate (on the (FC) form) \n', ...
        percent_c_d_FC, count_degenerate_FC );
    
    fprintf('\n');
    fprintf('\n');
    
    if nb_LCP_unsolved>0
        fprintf('the mean LCP error for unsolved LCP is: %.4e\n',mean(data_err_uns(:,1)));
    end
    
    fprintf('\n');
    count_degenerate = 0; count_degenerate_FC = 0;
    for i=1:nb_LCP_unsolved
        W = M_uns{i}*Z_uns{i} + Q_uns{i};
        supp_zer_sol = setdiff( 1:1:size(Z_uns{i},1) , ...
            union( find(abs(Z_uns{i})>1e-9),  find(abs(W)>1e-9) ) );
        if size(supp_zer_sol,2)>0
            count_degenerate = count_degenerate+1;
        end
        
        % for FC formulation (without \alpha-basis and only one
        % \beta-basis)
        dim_FC = size(Z_uns{i},1)/2; dim_nc = size(Z_uns{i},1)/4;
        E_tmp = union( find(abs(Z_uns{i})>1e-9),  find(abs(W)>1e-9) );
        E_tmp = E_tmp( E_tmp < 3*dim_nc+1 );
        E_tmp = E_tmp( E_tmp > dim_nc ); 
        rest_beta = setdiff( 1:1:dim_nc , ceil( (E_tmp-dim_nc) / 2 ) );
        rest_alpha = setdiff( 1:1:dim_nc , ...
            union( find( abs( Z_uns{i}(1:dim_nc) ) > 1e-9 ), ...
            find( abs( W(1:dim_nc) ) > 1e-9 ) ) );
        
        supp_zer_sol = union( rest_alpha , rest_beta );
        if size(supp_zer_sol,2)>0
            count_degenerate_FC = count_degenerate_FC+1;
        end
    end
    
    if nb_LCP_unsolved==0
        percent_c_d = 100;
    else
        percent_c_d     = count_degenerate      / nb_LCP_unsolved*100;
        percent_c_d_FC  = count_degenerate_FC   / nb_LCP_unsolved*100;
    end
    fprintf('%.4g%% (%d) of the unsolved LCP are degenerate (on the (APS) form) \n', ...
        percent_c_d, count_degenerate );
    fprintf('%.4g%% (%d) of the unsolved LCP are degenerate (on the (FC) form) \n', ...
        percent_c_d_FC, count_degenerate_FC );

%     unsolved_marker = 4;
%     lcp_error_status = 3;
%     member_name = strcat('/',GG_2,'/',Stats_solver);
%     data_unsol = h5read(file_name,member_name);
%     data_unsol = data_unsol'; % Why the reading is a transposed??
%     val_max_perturb = input('What is the maximum number of perturbation calls for solving? \n');
%     nb_algo = input('What is the maximum number of different solvers called? \n');
%     try
%         assert(all(data_unsol(:,6)==unsolved_marker),'The LCP is considere ''unsolved'' whereas the status from the lcp error is not ''unsolved''!')
%     catch
%         assert(all(data_unsol(:,1)==val_max_perturb),'The LCP is considere ''unsolved'' before reaching the maximum of perturbations');
%         assert(all(data_unsol(:,2)==nb_algo),'The LCP is considere ''unsolved'' before reaching the maximum of algorithms');
%     end
%     fprintf('The LCP is considered ''unsolved'' after %d perturbations of %d different solvers\n',val_max_perturb,nb_algo);
%     data_unsol = data_unsol(Idx_sol_alive_uns,:);
%     nb_diff_lcp_unsol = size(data_unsol,1);
%     
%     member_name = strcat('/',GG_1,'/',Stats_solver);
%     data_sol = h5read(file_name,member_name);
%     data_sol = data_sol'; % Why the reading is a transposed??
%     Algo_name_cpp = {'lexico','lemke'};
%     Algo_name_matlab = {'lemke','iterlemke','lexico'};
%     data_sol = data_sol(Idx_sol_alive,:);
%     nb_diff_lcp_sol = size(data_sol,1);
%     
%     per_algo = zeros(nb_algo,1);
%     for i=1:nb_algo
%         per_algo(i) = size(find(data_sol(:,2)==i),1)/...
%             nb_diff_lcp_sol*100;
%         fprintf('%.4g%% of successes are from the solver: %s\n',per_algo(i),Algo_name_cpp{i})
%     end
%     assert(abs(sum(per_algo)-100)<1e-10,'There exists algos do not taken into account')
%     
%     per_test = zeros(val_max_perturb+1,1);
%     for i=0:val_max_perturb
%         per_test(i+1) = size(find(data_sol(:,1)==i),1)/...
%             nb_diff_lcp_sol*100;
%         fprintf('%.4g%% of success after %d perturbations\n',per_test(i+1),i)
%     end
%     assert(abs(sum(per_test)-100)<1e-10,'There exists numbers of perturbations do not taken into account')
%     
%     for i=1:nb_algo
%         nb_ith_solver = sum( data_sol(data_sol(:,2)>=i,1) + 1 );
%         nb_other_solver = sum( data_sol(data_sol(:,2)<i,1) );
%         
%         per_algo = size(find(data_sol(:,2)==i),1);
%         fprintf('The solver: %s had a rate success of %.4g%%\n',Algo_name_cpp{i},...
%             per_algo/(nb_ith_solver+nb_other_solver)*100);
%     end
    
%     fprintf('\n');
%     if size(data_unsol,2)==6
%         per_status = zeros(lcp_error_status,1);
%         for i=1:lcp_error_status
%             per_status(i) = size(find(data_sol(:,6)==i),1)/nb_diff_lcp_sol*100;
%             fprintf('%.4g%% of solved LCP are with the status %d\n',per_status(i),i);
%         end
%         assert(abs(sum(per_status)-100)<1e-10,'There exists numbers of status do not taken into account');
%         
%         fprintf('\n');
%         
%         for i=1:nb_algo
%             success_algo = find(data_sol(:,2)==i);
%             per_status = zeros(lcp_error_status,1);
%             for j=1:lcp_error_status
%                 per_status(j) = size(find(data_sol(success_algo,6)==j),1)/...
%                     size(success_algo,1)*100;
%                 fprintf('%.4g%% of solved LCP with %s are with the status %d\n',...
%                     per_status(j),Algo_name_cpp{i},j);
%             end
%             assert(abs(sum(per_status)-100)<1e-10,'There exists numbers of status do not taken into account');
%         end
%         
%         fprintf('\n');
%         
%         for i=0:val_max_perturb
%             per_test = size(find(data_unsol(:,1)==i),1)/...
%                 nb_diff_lcp_unsol*100;
%             fprintf('The %dth perturbations is the best in %.4g%% of the unsolved cases\n',i,per_test);
%         end
%         
%         fprintf('\n');
%         
%         per_algo = zeros(nb_algo,1);
%         for i=1:nb_algo
%             per_algo(i) = size(find(data_unsol(:,2)==i),1)/...
%                 nb_diff_lcp_unsol*100;
%             fprintf('The solution from the solver %s is the best in %.4g%% of the unsolved cases\n',...
%                 Algo_name_cpp{i}, per_algo(i) );
%         end
%         assert(abs(sum(per_algo)-100)<1e-10,'There exists algos do not taken into account')
%         
%         fprintf('\n\n');
%         
%         fprintf('Which is the best way to perturb? (in term of success rate)\n')
%         diff_perturb = {'random','reduction'};
%         for i=1:size(diff_perturb,2)
%             Idx_perturb = find(data_sol(:,3)==i);
%             Idx_per_temp = data_sol(Idx_perturb,1)>=1;
%             Idx_per_used_sol = Idx_perturb(Idx_per_temp);
%             
%             Idx_perturb = find(data_unsol(:,3)==i);
%             Idx_per_temp = data_unsol(Idx_perturb,1)>=1;
%             Idx_per_used_unsol = Idx_perturb(Idx_per_temp);
%             
%             nb_total = size(Idx_per_used_sol,1) + size(Idx_per_used_unsol,1);
%             
%             success_rate = size(Idx_per_used_sol,1)/nb_total*100;
%             mean_nb_call = mean(data_sol(Idx_per_used_sol,1));
%             mean_qual_pert = mean(data_sol(Idx_per_used_sol,6));
%                         
%             fprintf('The %s perturbation is called %d times, has a success rate of %.4g%%, needs %.4g call in average, with a quality of %.4g\n',...
%                 diff_perturb{i}, nb_total, success_rate, mean_nb_call,...
%                 mean_qual_pert );
%             
%             fprintf('\n');    
%         end
%     end
%     
%     member_name = strcat('/',GG_2,'/',Stats_err);
%     data_err_uns = h5read(file_name,member_name);
%     member_name = strcat('/',GG_1,'/',Stats_err);
%     data_err = h5read(file_name,member_name);
%     
%     data_err_uns = data_err_uns(Idx_sol_alive_uns);    
%     data_err = data_err(Idx_sol_alive);
    
    %% saving data:
    fprintf('Saving phase...\n');
    mat_save_name = strcat('mat_',suffix{1});
    if strcmp(ans_rem_sim_LCP,'y')
        save(mat_save_name,'M','Q','Z','Ep','M_uns','Q_uns','Z_uns','Ep_uns',...
        'is_same_LCP','Idx_diff_uns_sol','Idx_sol_alive','Idx_sol_alive_uns',...
        'list_same_lcp_sol','list_same_lcp_unsol','data_sol','data_unsol',...
        'data_err','data_err_uns','Contact_Graph');
    else
        save(mat_save_name,'M','Q','Z','M_uns','Q_uns','Z_uns',...
            'data_err','data_err_uns','Contact_Graph');
    end
    
    %% SOL with Lemke, Iterlemke and lexicolemke solvers:
    fprintf('Resolution with solvers from Matlab...\n');
%     SOL_lemke_uns = cell(size(M_uns,1),3);
%     SOL_IterLemke_uns = cell(size(M_uns,1),4);
%     SOL_lexico_uns = cell(size(M_uns,1),3);
%     
%     SOL_lemke = cell(size(M,1),3);
%     SOL_IterLemke = cell(size(M,1),4);
%     SOL_lexico = cell(size(M,1),3);
    
    SOL_Geo_lexi = cell(size(M,1),6);
    SOL_Geo_lexi_uns = cell(size(M_uns,1),6);
    
%     t = cputime;
    count_solved = 0; % err_lexico = -1*ones(size(M,1),1);
%     err_lcp_lexi = zeros(size(M,1),1); iter_lcp = err_lexico; cone_chgt = iter_lcp;
    % use_lexico_order = err_lcp_lexi;
    list_uns_sol_matlab = []; list_sol_SR_matlab = []; list_lexico_cmp = [];
    for i=1:size(M,1)
        if mod(i,500)==0; fprintf('%.4g%%\n',i/size(M,1)*100); end
        
%         %% scaling (first attempt):
%         dim = size(M{i},1);
%         DM = M{i}(1:3*dim/4,1:3*dim/4);
%         DM_wo_zero = DM(:);
%         scal_coef = 1/max(abs(DM_wo_zero));
%         DM_tilde = scal_coef*DM;
%         Q_tilde = scal_coef*Q{i};
%         M_tilde = [DM_tilde M{i}(1:3*dim/4,3*dim/4+1:end); ...
%             M{i}(3*dim/4+1:end,:)];
%         %% end of scaling
        M_tilde = M{i}; Q_tilde = Q{i};
        
        [z_lex,info] = Geo_lexico_Lemke(M_tilde,Q_tilde);
        W_lex = M{i}*z_lex+Q{i};
        err_lcp = sum(abs(z_lex.*W_lex)) + sum(abs(z_lex(z_lex<0))) + ...
            sum(abs(W_lex(W_lex<0)));assert(abs(err_lcp-info.err)<1e-7,'different error computation?!!?');
        if err_lcp<=1e-7
            count_solved = count_solved+1;
        else
            list_uns_sol_matlab = [list_uns_sol_matlab ; i];
        end
        
        SOL_Geo_lexi{i,1} = z_lex;
        SOL_Geo_lexi{i,2} = info.err_idx; SOL_Geo_lexi{i,3} = info.tot_iter;
        SOL_Geo_lexi{i,4} = info.err;
        SOL_Geo_lexi{i,5} = info.nb_adj_cone_covered;
        SOL_Geo_lexi{i,6} = info.tol;
        
        if info.nb_adj_cone_covered>0
            list_sol_SR_matlab = [list_sol_SR_matlab ; i];
        end
        
        if info.used_lexico>0
            list_lexico_cmp = [list_lexico_cmp ; i];
        end
    end
    
    fprintf('It stay only: %d LCP no solved amongst solved one.\n',size(M,1)-count_solved);
    
    count_solved = 0; list_uns_unsol_matlab = [];
    list_unsol_SR_matlab = []; list_lexico_cmp_uns = [];
    for i=1:size(M_uns,1)
        if mod(i,500)==0; fprintf('%.4g%%\n',i/size(M_uns,1)*100); end
        
%         %% scaling (first attempt):
%         dim = size(M_uns{i},1);
%         DM = M_uns{i}(1:3*dim/4,1:3*dim/4);
%         DM_wo_zero = DM(:);
%         scal_coef = 1/max(abs(DM_wo_zero));
%         DM_tilde = scal_coef*DM;
%         Q_tilde_uns = scal_coef*Q_uns{i};
%         M_tilde_uns = [DM_tilde M_uns{i}(1:3*dim/4,3*dim/4+1:end); ...
%             M_uns{i}(3*dim/4+1:end,:)];
%         %% end of scaling
        M_tilde_uns = M_uns{i}; Q_tilde_uns = Q_uns{i};
 
        [z_lex,info] = Geo_lexico_Lemke(M_tilde_uns,Q_tilde_uns);
        W_lex = M_uns{i}*z_lex+Q_uns{i};

        err_lcp = sum(abs(z_lex.*W_lex)) + sum(abs(z_lex(z_lex<0))) + ...
            sum(abs(W_lex(W_lex<0)));assert(abs(err_lcp-info.err)<1e-7,'different error computation?!!?');

        if err_lcp<=1e-7
            count_solved = count_solved+1;
        else
            list_uns_unsol_matlab = [list_uns_unsol_matlab ; i];
        end

        SOL_Geo_lexi_uns{i,1} = z_lex;
        SOL_Geo_lexi_uns{i,2} = info.err_idx; SOL_Geo_lexi_uns{i,3} = info.tot_iter;
        SOL_Geo_lexi_uns{i,4} = info.err;
        SOL_Geo_lexi_uns{i,5} = info.nb_adj_cone_covered;
        SOL_Geo_lexi_uns{i,6} = info.tol;

        if info.nb_adj_cone_covered>0
            list_unsol_SR_matlab = [list_unsol_SR_matlab ; i];
        end
        if info.used_lexico>0
            list_lexico_cmp_uns = [list_lexico_cmp_uns ; i];
        end
    end
    
    fprintf('It stay only: %d LCP no solved amongst unsolved one.\n',size(M_uns,1)-count_solved);
   
    fprintf('%d LCP solver (solved one) ends on the secondary ray\n', ...
        size(list_sol_SR_matlab,1));
    fprintf('%d LCP solver (unsolved one) ends on the secondary ray\n', ...
        size(list_unsol_SR_matlab,1));
    
    fprintf('\n\n');
    fprintf('%d LCP solver (solved one) used the lexicographical ordering\n', ...
        size(list_lexico_cmp,1));
    fprintf('%d LCP solver (unsolved one) used the lexicographical ordering\n', ...
        size(list_lexico_cmp_uns,1));
    
    var_tol = zeros(size(M,1),2); var_tol_uns = zeros(size(M_uns,1),2);
    LCP_err = zeros(size(M,1),1); LCP_err_uns = zeros(size(M_uns,1),1);
    for i=1:size(M,1)
        var_tol(i,:) = SOL_Geo_lexi{i,6}; LCP_err(i) = SOL_Geo_lexi{i,4};
    end
    for i=1:size(M_uns,1)
        var_tol_uns(i,:) = SOL_Geo_lexi_uns{i,6}; LCP_err_uns(i) = SOL_Geo_lexi_uns{i,4};
    end
    fprintf('\n\n');
    fprintf('%d Relative velocities contains a null coefficient (amongst solved one)\n', ...
        size( find( var_tol(:,1)==1 ), 1) );
    fprintf('%d times the tolerance has been modified (amongst solved one)\n', ...
        size( find( var_tol(:,2)==1 ), 1) );
    fprintf('\n');
    fprintf('%d Relative velocities contains a null coefficient (amongst unsolved one)\n', ...
        size( find( var_tol_uns(:,1)==1 ), 1) );
    fprintf('%d times the tolerance has been modified (amongst unsolved one)\n', ...
        size( find( var_tol_uns(:,2)==1 ), 1) );
    
    set_of_zeros_RVel = find( var_tol(:,1)==1 );
    set_of_mod_tol = find( var_tol(:,2)==1 );
    set_of_zeros_RVel_uns = find( var_tol_uns(:,1)==1 );
    set_of_mod_tol_uns = find( var_tol_uns(:,2)==1 );
    temp_Idx_RV = setdiff(list_uns_sol_matlab,set_of_zeros_RVel);
    temp_Idx_RV_uns = setdiff(list_uns_unsol_matlab,set_of_zeros_RVel_uns);
    temp_Idx_MT = setdiff(list_uns_sol_matlab,set_of_mod_tol);
    temp_Idx_MT_uns = setdiff(list_uns_unsol_matlab,set_of_mod_tol_uns);
    
    if size(set_of_zeros_RVel,1)==0
        per_RelV_uns = 0;
    else
        per_RelV_uns = abs(size(temp_Idx_RV,1) - size(list_uns_sol_matlab,1))/...
            size(set_of_zeros_RVel,1)*100;
    end
    if size(set_of_zeros_RVel_uns,1)==0
        per_RelV_uns_uns = 0;
    else
        per_RelV_uns_uns = abs(size(temp_Idx_RV_uns,1) - size(list_uns_unsol_matlab,1))/...
            size(set_of_zeros_RVel_uns,1)*100;
    end
    
    if size(set_of_mod_tol,1)==0
        per_MT_uns = 0;
    else
        per_MT_uns = abs(size(temp_Idx_MT,1) - size(list_uns_sol_matlab,1))/...
            size(set_of_mod_tol,1)*100;
    end
    if size(set_of_mod_tol_uns,1)==0
        per_MT_uns_uns = 0;
    else
        per_MT_uns_uns = abs(size(temp_Idx_MT_uns,1) - size(list_uns_unsol_matlab,1))/...
            size(set_of_mod_tol_uns,1)*100;
    end
    
    fprintf('\n\n');
    fprintf('%.4g%% Relative velocities containing a null coefficient (amongst solved one) are unsolved\n', ...
         per_RelV_uns);
    fprintf('%.4g%% Relative velocities containing a null coefficient (amongst unsolved one) are unsolved\n', ...
         per_RelV_uns_uns);
    fprintf('\n');
    fprintf('%.4g%% Modified tolerance (amongst solved one) are unsolved\n', ...
         per_MT_uns);
    fprintf('%.4g%% Modified tolerance (amongst unsolved one) are unsolved\n', ...
         per_MT_uns_uns);
%     e_lex = cputime - t;

%     mean(err_lcp_lexi)
%     size(find(cone_chgt~=0))
%     mean(iter_lcp)
%     max(iter_lcp)
%     size(find(err_lexico==0))


    
    % LCP error: 3 part : Energy part (z^Tw), velocity (relative) part (w) and impulse
    % part (z).
    % We could set limits compared to some physical properties form sea
    % ice: 
    % 1/ maximal velocity (1m/s) => abs(w) < 2m/s (relative velocity)
    %    acceptable error from abs(w) is 1e-5.
    % 2/ nothing on energy => transform to velocity : V = sqrt( z^Tw / M)
    %    where M is the mass of the ice floe. Thus, V < 1m/s, therefore
    %    abs(z^T .* w) < ( 5e-6 )^2 * M 
    %    From the capability of the model, I set to min(M)=1e5 kg and
    %    max(M) = 1e11kg (To this day: 1 -> 100 for the size, 1 -> 1000 is
    %    not reached yet).
    %    Thus acceptable error from abs(z^T .* w) is 2.5e-6.
    %    Let us imagine an application for interaction ice/structures. We
    %    need to reduce size until 1 m (even 0.5 m) => 5e2 kg.
    %    In this case, velocity could be greater and we could increase the
    %    acceptable error in velocity to 5e-5 => leading to acceptable
    %    error from abs(z^T .* w) set to 1.25e-6.
    % 3/ part from contact impulse (real values unknown). Acceptable error
    %    taken into account in 1/ and 2/.
    %
    % Two physical limits are important: 1/ for non interpenetration and 2/
    % for no increase of energy due to contact.
    % These two limits could be taken into account with acceptable error on
    % the first part of abs(w) (corresponding to normal relative
    % velocities) and on the second point (error on abs(z^T .* w)).
    % 
    % Let us setted up the first to 1e-6 and the second to 1e-7.
    % Observing the efficiency of the solver to give solution satisfying
    % these physical limits.
    tol_err_p1 = 1e-6;
    tol_err_p2 = 1e-7;
%     fprintf('The acceptable LCP error are: %.1e (corresponding to the relative normal velocity) and %.1e (corresponding to the increase of kinetic energy during the contact).\n',...
%         tol_err_p1 , tol_err_p2 );
%     
%     Mat_Err = []; Idx_max_Err = []; Mass = -1;
%     for i=1:size(M_uns,1)
%         A = M_uns{i}; q = Q_uns{i};
%         nb_c = size(A,1)/4;
%         
%         % LEMKE:
%         [z_lem, err, iter] = lemke(A,q);
%         SOL_lemke_uns{i,1} = z_lem;
%         W_lem = A*z_lem+q;
%         err_lcp = sum(abs(z_lem.*W_lem)) + sum(abs(z_lem(z_lem<0))) + ...
%             sum(abs(W_lem(W_lem<0)));
%         SOL_lemke_uns{i,2} = err_lcp; SOL_lemke_uns{i,3} = iter;
%         SOL_lemke_uns{i,4} = err;
%         
%         % IterLEMKE:
%         [z_I, ~, info] = IterLemke(A,q);
%         SOL_IterLemke_uns{i,1} = z_I;
%         W_I = A*z_I+q;
%         SOL_IterLemke_uns{i,2} = info.sol_err; SOL_IterLemke_uns{i,3} = info.lmk_iter;
%         SOL_IterLemke_uns{i,4} = info.error_code;
%         
%         % LexicoLEMKE:
%         [z_lex, ~, iter] = lexicolemkeplusopt(A,q,10000);
%         SOL_lexico_uns{i,1} = z_lex;
%         W_lex = A*z_lex+q;
%         err_lcp = sum(abs(z_lex.*W_lex)) + sum(abs(z_lex(z_lex<0))) + ...
%             sum(abs(W_lex(W_lex<0)));
%         SOL_lexico_uns{i,2} = err_lcp; SOL_lexico_uns{i,3} = iter;
%         
%         % Mat_Err:
%         [Mat_Err,Idx_max_Err] = Calc_Mat_Err(0,Mat_Err,Idx_max_Err,z_lem,W_lem,...
%             nb_c,Mass,[tol_err_p1,tol_err_p2],i,1);
%         [Mat_Err,Idx_max_Err] = Calc_Mat_Err(0,Mat_Err,Idx_max_Err,z_I,W_I,...
%             nb_c,Mass,[tol_err_p1,tol_err_p2],i,2);
%         [Mat_Err,Idx_max_Err] = Calc_Mat_Err(0,Mat_Err,Idx_max_Err,z_lex,W_lex,...
%             nb_c,Mass,[tol_err_p1,tol_err_p2],i,3);
%     end
%     
%     for i=1:size(M,1)
%         if mod(i,500)==0; fprintf('%.4g%%\n',i/size(M,1)*100); end
%         A = M{i}; q = Q{i};
%         nb_c = size(A,1)/4;
%                 
%         % LEMKE:
%         [z_lem, err, iter] = lemke(A,q);
%         SOL_lemke{i,1} = z_lem;
%         W_lem = A*z_lem+q;
%         err_lcp = sum(abs(z_lem.*W_lem)) + sum(abs(z_lem(z_lem<0))) + ...
%             sum(abs(W_lem(W_lem<0)));
%         SOL_lemke{i,2} = err_lcp; SOL_lemke{i,3} = iter;
%         SOL_lemke{i,4} = err;
%         
%         % IterLEMKE:
%         [z_I, ~, info] = IterLemke(A,q);
%         SOL_IterLemke{i,1} = z_I;
%         W_I = A*z_I+q;
%         SOL_IterLemke{i,2} = info.sol_err; SOL_IterLemke{i,3} = info.lmk_iter;
%         SOL_IterLemke{i,4} = info.error_code;
%         
%         % LexicoLEMKE:
%         [z_lex, ~, iter] = lexicolemkeplusopt(A,q,10000);
%         SOL_lexico{i,1} = z_lex;
%         W_lex = A*z_lex+q;
%         err_lcp = sum(abs(z_lex.*W_lex)) + sum(abs(z_lex(z_lex<0))) + ...
%             sum(abs(W_lex(W_lex<0)));
%         SOL_lexico{i,2} = err_lcp; SOL_lexico{i,3} = iter;
%         
%         % Mat_Err:
%         [Mat_Err,Idx_max_Err] = Calc_Mat_Err(1,Mat_Err,Idx_max_Err,z_lem,W_lem,...
%             nb_c,Mass,[tol_err_p1,tol_err_p2],i,1);
%         [Mat_Err,Idx_max_Err] = Calc_Mat_Err(1,Mat_Err,Idx_max_Err,z_I,W_I,...
%             nb_c,Mass,[tol_err_p1,tol_err_p2],i,2);
%         [Mat_Err,Idx_max_Err] = Calc_Mat_Err(1,Mat_Err,Idx_max_Err,z_lex,W_lex,...
%             nb_c,Mass,[tol_err_p1,tol_err_p2],i,3);
%     end
%     
%     % Success rate for Matlab solvers:
%     if nb_LCP_solved~=0 || nb_LCP_unsolved~=0
%         for i=1:3
%             SR = size(find(Mat_Err(:,18)==i),1);
%             fprintf('The success rate of %s is %.4g%%.\n',Algo_name_matlab{i},...
%             100-SR/(nb_LCP_solved+nb_LCP_unsolved)*100);
%         end
%         
%         List_contact_unsolved = [];
%         Ens_matlab = Mat_Err(:,2);
%         Ens_matlab = unique(Ens_matlab);
%         for i=1:size(Ens_matlab,1)
%             Ens_temp = Mat_Err(:,2)==Ens_matlab(i);
%             idx_temp = setdiff([1 2 3]',Mat_Err(Ens_temp,18));
%             if size(idx_temp,1)==0
%                 List_contact_unsolved = [List_contact_unsolved;Ens_matlab(i)];
%             end
%         end
%         
%         fprintf('There are %.4g%% (i.e.: %d) of unsolved LCP by any Matlab solvers\n',...
%             size(List_contact_unsolved,1)/(nb_LCP_solved+nb_LCP_unsolved)*100,...
%             size(List_contact_unsolved,1));
%     end
%     
    %% Comparison with solution and error from C++ implementation:
    fprintf('Checking phase for solutions from C++ solvers...\n');
    % checking LCP error:
%     diff_err = zeros(nb_LCP_unsolved,1); 
%     Mat_Err_cpp = []; Idx_max_Err_cpp = []; Mass = -1;
%     for i=1:nb_LCP_unsolved
%         z = Z_uns{i}; W = M_uns{i}*z + Q_uns{i};
%         [~,~,~,e] = lcp_error(z,W); 
%         nb_c = size(z,1)/4;
%         
%         [Mat_Err_cpp,Idx_max_Err_cpp] = Calc_Mat_Err(0,Mat_Err_cpp,Idx_max_Err_cpp,...
%             z,W,nb_c,Mass,[tol_err_p1,tol_err_p2],i,0);
%         
%         if data_err_uns(i)>1
%             diff_err(i) = (e-data_err_uns(i))/data_err_uns(i);
%         else
%             diff_err(i) = e-data_err_uns(i);
%         end
%     end
%     if size(diff_err,1)>0
%         assert(max(abs(diff_err))<1e-7,'Bigre!!!')
%     end
%     
%     diff_err = zeros(nb_LCP_solved,1);
%     for i=1:nb_LCP_solved
%         if mod(i,500)==0; fprintf('%.4g%%\n',i/size(M,1)*100); end
%         z = Z{i}; W = M{i}*z + Q{i};
%         [~,~,~,e] = lcp_error(z,W); 
%         nb_c = size(z,1)/4;
%         
%         [Mat_Err_cpp,Idx_max_Err_cpp] = Calc_Mat_Err(1,Mat_Err_cpp,Idx_max_Err_cpp,...
%             z,W,nb_c,Mass,[tol_err_p1,tol_err_p2],i,0);
%         
%         if data_err(i)>1
%             diff_err(i) = (e-data_err(i))/data_err(i);
%         else
%             diff_err(i) = e-data_err(i);
%         end
%     end
%     if size(diff_err,1)>0
%         assert(max(abs(diff_err))<1e-7,'Bigre!!!')
%     end
    
    fprintf('\n');
    % checking solution (z):
%     diff_Z_uns = zeros(nb_LCP_unsolved,3); diff_Z = zeros(nb_LCP_solved,3);
    diff_Z_uns = zeros(nb_LCP_unsolved,1); diff_Z = zeros(nb_LCP_solved,1);
    for i=1:nb_LCP_unsolved
%         z_lem = SOL_lemke_uns{i,1} - Z_uns{i};
%         z_I = SOL_IterLemke_uns{i,1} - Z_uns{i};
%         z_lex = SOL_lexico_uns{i,1} - Z_uns{i};
%         diff_Z_uns(i,:) = [norm(z_lem) norm(z_I) norm(z_lex)];
        
        z_geo = SOL_Geo_lexi_uns{i,1} - Z_uns{i};
        diff_Z_uns(i) = norm(z_geo);
    end
    
    for i=1:nb_LCP_solved
        if mod(i,500)==0; fprintf('%.4g%%\n',i/size(M,1)*100); end
%         z_lem = SOL_lemke{i,1} - Z{i};
%         z_I = SOL_IterLemke{i,1} - Z{i};
%         z_lex = SOL_lexico{i,1} - Z{i};
%         diff_Z(i,:) = [norm(z_lem) norm(z_I) norm(z_lex)];
        
        z_geo = SOL_Geo_lexi{i,1} - Z{i};
        diff_Z(i) = norm(z_geo);
    end
    
    if nb_LCP_solved>0
        fprintf('%.4g%% (%d) of the solutions from solved cases are similar to Geo_info in Matlab \n',...
            size(find(diff_Z==0),1)/nb_LCP_solved*100, size(find(diff_Z==0),1));
        if size(find(diff_Z~=0),1)>0
            fprintf('The mean difference between solutions not similar is: %.4e \n',...
                mean(diff_Z(diff_Z~=0))/size(find(diff_Z~=0),1)*1);
        end
    end
    if nb_LCP_unsolved>0
        fprintf('%.4g%% (%d) of the solutions from unsolved cases are similar to Geo_info in Matlab\n',...
            size(find(diff_Z_uns==0),1)/nb_LCP_unsolved*100, size(find(diff_Z_uns==0),1));
        if size(find(diff_Z~=0),1)>0
            fprintf('The mean difference between solutions not similar is: %.4e \n',...
                mean(diff_Z(diff_Z~=0))/size(find(diff_Z_uns~=0),1)*1);
        end
    end

    fprintf('\n\n')
    if size(list_uns_sol_matlab,1)>0
        miss_matlab = setdiff(list_uns_sol_matlab,find(data_err(:,3)>0));
        if size(miss_matlab,1)==0 || size(miss_matlab,2)==0
            s_m_m = 0;
        else
            s_m_m = size(miss_matlab,1);
        end
        if size(find(data_err(:,3)>0),1)
            fprintf('%.4g%% (%d) of the unsolved LCP from matlab routines are solved in c++ by perturbations \n',...
                ( size(find(data_err(:,3)>0),1) - s_m_m )/...
                size(find(data_err(:,3)>0),1)*100, ...
                ( size(find(data_err(:,3)>0),1) - s_m_m ) );
        end
    end
    
    if size(list_sol_SR_matlab,1)>0
        miss_matlab = setdiff(list_sol_SR_matlab,find(data_err(:,4)>0));
        if size(miss_matlab,1)==0 || size(miss_matlab,2)==0
            s_m_m = 0;
        else
            s_m_m = size(miss_matlab,1);
        end
        if size(find(data_err(:,4)>0),1)
            fprintf('%.4g%% (%d) of the solved LCP from matlab routines with SR are also solved in c++ with SR \n',...
                ( size(find(data_err(:,4)>0),1) - s_m_m )/...
                size(find(data_err(:,4)>0),1)*100, ...
                ( size(find(data_err(:,4)>0),1) - s_m_m ) );
        end
    end
    
    if size(list_lexico_cmp,1)>0
        miss_matlab = setdiff(list_lexico_cmp ,find(data_err(:,6)>0));
        if size(miss_matlab,1)==0 || size(miss_matlab,2)==0
            s_m_m = 0;
        else
            s_m_m = size(miss_matlab,1);
        end
        if size(find(data_err(:,6)>0),1)
            fprintf('%.4g%% (%d) of the solved LCP from matlab routines with lexico ordering are also solved with lexico ordering \n',...
                ( size(find(data_err(:,6)>0),1) - s_m_m )/...
                size(find(data_err(:,6)>0),1)*100, ...
                ( size(find(data_err(:,6)>0),1) - s_m_m ) );
        end
    end
    
%     int_temp = find(data_sol(:,2)==2);
%     for i=1:3
%         fprintf('%.4g%% of the solutions from unsolved cases are similar to %s in Matlab\n',...
%             size(find(diff_Z_uns(:,i)==0),1)/nb_LCP_unsolved*100,Algo_name_matlab{i});
%         fprintf('%.4g%% of the solutions from solved cases are similar to %s in Matlab\n',...
%             size(find(diff_Z(:,i)==0),1)/nb_LCP_solved*100,Algo_name_matlab{i});
%         fprintf('%.4g%% of the solutions from lemke in C++ are similar to %s in Matlab\n',...
%             size(find(diff_Z(int_temp,i)==0),1)/size(int_temp,1)*100,Algo_name_matlab{i});
%         fprintf('%.4g%% of the solutions from lemke in C++ are similar to %s in Matlab more or less 1e-10\n',...
%             size(find(diff_Z(int_temp,i)<1e-10),1)/size(int_temp,1)*100,Algo_name_matlab{i});
%         fprintf('\n\n');
%     end
    
    % Analysing of the source of errors and comparison with Matlab solvers:
    if nb_LCP_unsolved>0
    rest_100 = mod(data_err_uns(:,end-1),100);
    rest_3 = mod(rest_100,20);
    rest_20 = rest_100-rest_3;
    
    fprintf('%.4g%% of the unsolved LCP are due to relative normale velocity that may cause an interpenetration\n',...
        size(find(rest_3==3),1)/nb_LCP_unsolved*100);
    fprintf('%.4g%% of the unsolved LCP are due to an increase of kinetic energy\n',...
        size(find(rest_20==20),1)/nb_LCP_unsolved*100);
    fprintf('%.4g%% of the unsolved LCP are due to impulse error above the limit\n',...
        size(find(data_err_uns(:,end)>=100),1)/nb_LCP_unsolved*100);
    fprintf('\n\n');
    
    source_part1 = mod(Mat_Err_cpp(:,16),10);
    source_part2 = (Mat_Err_cpp(:,16) - source_part1)/10;
    nc = Mat_Err_cpp(:,3);
    source_part = 3*ones(size(source_part1,1),2);
    source_part(:,2) = source_part1;
    source_part(source_part2<=3*nc,1) = 2;
    source_part(source_part2<=nc,1) = 1;

    source_text = {'\lambda^T J^T W^+','\beta^T D^T W^+',...
        '\alpha^T \left(\mu\lambda-H^T \beta\right)','\lambda','\beta',...
        '\alpha','J^T W^+','D^T W^+','contact forces outside of the cone'};
    
    for i=1:3
        A = source_part(:,1)==i;
        for j=1:3
            B = source_part(:,2)==j;
            C = A.*B;
            k = (i-1)*3+j;
            
            fprintf('The largest part of the LCP error comes from %s in %.4g%% of cases\n',...
                source_text{k},size(find(C==1),1)/size(Mat_Err_cpp,1)*100)
        end
    end
    fprintf('\n\n');
    
    % Hypothesis: z from Lemke_eigen in C++ in the unsolved case is in the
    % form: \lambda, \beta << \alpha!
    count = 0; list_hyp_false = [];
    for i=1:nb_LCP_unsolved
        z = Z_uns{i}; nc = size(z,1)/4;
        z_alpha = max(abs(z(3*nc+1:end)));
        z_other = max(abs(z(1:3*nc)));
        if z_alpha~=0
            if z_other/z_alpha<1e-7
                count = count+1;
            else
                list_hyp_false = [list_hyp_false;i];
            end
        end
    end
    fprintf('%.4g%% of the unsolved cases are due to unrealistic coefficient alpha\n',...
        count/nb_LCP_unsolved*100);
    end
    fprintf('\n\n');
    
    % Checking solutions from Matlab solvers satisfy physical conditions
    % and LCP error (Relative normal velocity > 0 and no increase of 
    % kinetic energy):
%     List_contact_unsolved_cpp = [];
%     for i=1:size(List_contact_unsolved,1)
%         idx_temp = find(Mat_Err_cpp(:,2)==List_contact_unsolved(i));
%         if size(idx_temp,1)~=0
%             List_contact_unsolved_cpp = [List_contact_unsolved_cpp;...
%                 List_contact_unsolved(i)];
%         end
%     end
%     
%     fprintf('There are %.4g%% (i.e.: %d) of unsolved LCP by any solvers from Matlab and C++\n',...
%         size(List_contact_unsolved_cpp,1)/(nb_LCP_solved+nb_LCP_unsolved)*100,...
%             size(List_contact_unsolved_cpp,1));
    
    SOL_save_name = strcat('SOL_',suffix{1});
    save(SOL_save_name,'SOL_Geo_lexi','SOL_Geo_lexi_uns',...
        'Mat_Err_cpp','Idx_max_Err_cpp','diff_Z','diff_Z_uns',...
        'list_uns_sol_matlab','list_sol_SR_matlab',...
        'list_sol_SR_matlab','list_unsol_SR_matlab',...
        'list_lexico_cmp','list_lexico_cmp_uns');
        
%     save(SOL_save_name,'SOL_lemke','SOL_IterLemke','SOL_lexico',...
%         'SOL_lemke_uns','SOL_IterLemke_uns','SOL_lexico_uns','Mat_Err',...
%         'Idx_max_Err','Mat_Err_cpp','Idx_max_Err_cpp',...
%         'List_contact_unsolved','List_contact_unsolved_cpp');
    
    bool_save_FC_form = input('Do you want save ''M'' and ''q'' from LCP(M,q) on the convex-analysis form?\n 0 for No, 1 for Yes\n');
    
    %% Saving on the FC form (removing facetisation of Coulomb's Cone + 
    % removing the last part due to Lagrangian coefficient (\alpha):
    if bool_save_FC_form
        M_tilde = cell(nb_LCP_solved,1); Q_tilde = M_tilde;
        M_uns_tilde = cell(nb_LCP_unsolved,1); Q_uns_tilde = M_uns_tilde;
        
        for i=1:nb_LCP_solved
            nbc = size(M{i},1)/4;
            assert(floor(nbc)==nbc,'BIGRE!!!')
            s_del = 3*nbc;
            ind_rl = [1:nbc nbc+1:2:s_del];
            
            M_tilde{i} = M{i}(ind_rl,ind_rl); Q_tilde{i} = Q{i}(ind_rl);
        end
        for i=1:nb_LCP_unsolved
            nbc = size(M_uns{i},1)/4;
            assert(floor(nbc)==nbc,'BIGRE!!!')
            s_del = 3*nbc;
            ind_rl = [1:nbc nbc+1:2:s_del];
            
            M_uns_tilde{i} = M_uns{i}(ind_rl,ind_rl); Q_uns_tilde{i} = Q_uns{i}(ind_rl);
        end
    end
    
    %% Testing with reduction of LCP in the form M + \alpha I, with:
    % I: the identity matrix
    % \alpha: a coefficient lesser than 1 
    % The property: the LCP(M,q) is reductible is not valid here because M
    % is only co-positive, however, trying to compare the efficiency of
    % this method against a random perturbation:
    Algo_name_matlab = {'LexicoLemke_MR'};
    s_algo_matlab = size(Algo_name_matlab,2);
    coef_perturb = 1e-9; ite_max_perturb = 5;
    red_term = coef_perturb*ones(ite_max_perturb,1);
    perturb1 = coef_perturb; perturb2 = coef_perturb;
    
    fprintf('The total coefficient for the reduction is %.2e\n',sum(red_term));
    fprintf('The perturbation for random1 is %.2e, and for random2 is %.2e\n',...
        perturb1,perturb2);
    fprintf('The max of perturbation is %d\n',ite_max_perturb);
    fprintf('\n');
    max_count = size(red_term,2);
    
    tol_lcp = min(tol_err_p1,tol_err_p2);
    ite_lcp = cell(s_algo_matlab,1);
    Idx_solver = [1:1:nb_LCP_solved+nb_LCP_unsolved]';
    for i=1:s_algo_matlab
%         Idx_solver = find(Mat_Err(:,end)==i);
        ite_lcp{i} = zeros(size(Idx_solver,1),9);
        % ite_lcp{i}(j,:) = [nb ite to find with test1 | best_err | at count | ...
        % nb ite to find with test2 | best_err | at count | ...
        % nb ite to find with test3 | best_err | at count];
        
        min_coef = zeros(size(Idx_solver));
        for j=1:size(Idx_solver,1)
            if mod(j,50)==0; fprintf('%.4g%%\n',j/size(Idx_solver,1)*100); end
%             if Mat_Err(Idx_solver(j),1)==0 % from unsolved
            if j<= nb_LCP_solved % from solved
                M_orig = M{Idx_solver(j)};
                q_orig = Q{Idx_solver(j)};
            else % from unsolved //from solved
%                 M_orig = M{Mat_Err(Idx_solver(j),2)};
%                 q_orig = Q{Mat_Err(Idx_solver(j),2)};
                M_orig = M_uns{Idx_solver(j)-nb_LCP_solved};
                q_orig = Q_uns{Idx_solver(j)-nb_LCP_solved};
            end
            
            A = M_orig; best_err = Inf; best_count = 0;
            alpha = perturb1; solved = 0; count = 1;
            % first random
            while ~solved && count <= max_count
                A = random_perturb(alpha,A);
                    
                switch i
                    case 1 % LEMKE:
                        [z,~] = Geo_lexico_Lemke(A,q_orig); 
%                         [z, ~, ~] = lemke(A,q_orig);
                    case 2 % IterLEMKE:
                        [z, ~, ~] = IterLemke(A,q_orig);
                    case 3 % LexicoLEMKE:
                        [z, ~, ~] = lexicolemkeplusopt(A,q_orig,10000);
                end
                W = M_orig*z+q_orig;
                err_lcp = sum( abs(z.*W) ) + sum( abs( z( z<0 ) ) ) + ...
                    sum( abs( W( W<0 ) ) );
                
                if err_lcp <= tol_lcp
                    solved = 1; ite_lcp{i}(j,1) = count;
                end
                if err_lcp < best_err
                    best_err = err_lcp; best_count = count;
                end
                count = count+1;
            end
            ite_lcp{i}(j,2:3) = [best_err, best_count];
            
            A = M_orig; best_err = Inf; best_count = 0;
            alpha = perturb2; solved = 0; count = 1;
            % second random
            while ~solved && count <= max_count
                A = random_perturb2(alpha,A);
                  
                switch i
                    case 1 % LEMKE:
                        [z,~] = Geo_lexico_Lemke(A,q_orig);
%                         [z, ~, ~] = lemke(A,q_orig);
                    case 2 % IterLEMKE:
                        [z, ~, ~] = IterLemke(A,q_orig);
                    case 3 % LexicoLEMKE:
                        [z, ~, ~] = lexicolemkeplusopt(A,q_orig,10000);
                end
                W = M_orig*z+q_orig;
                err_lcp = sum( abs(z.*W) ) + sum( abs( z( z<0 ) ) ) + ...
                    sum( abs( W( W<0 ) ) );
                
                if err_lcp <= tol_lcp
                    solved = 1; ite_lcp{i}(j,4) = count;
                end
                if err_lcp < best_err
                    best_err = err_lcp; best_count = count;
                end
                count = count+1;
            end 
            ite_lcp{i}(j,5:6) = [best_err, best_count];
            
            A = M_orig; best_err = Inf; best_count = 0;
            solved = 0; count = 1;
            % third test:
            while ~solved && count <= max_count
                alpha = red_term(count);
                [A,min_coef(j)] = reduction_I(alpha,A);

                switch i
                    case 1 % LEMKE:
                        [z,~] = Geo_lexico_Lemke(A,q_orig);
%                         [z, ~, ~] = lemke(A,q_orig);
                    case 2 % IterLEMKE:
                        [z, ~, ~] = IterLemke(A,q_orig);
                    case 3 % LexicoLEMKE:
                        [z, ~, ~] = lexicolemkeplusopt(A,q_orig,10000);
                end
                W = M_orig*z+q_orig;
                err_lcp = sum( abs(z.*W) ) + sum( abs( z( z<0 ) ) ) + ...
                    sum( abs( W( W<0 ) ) );
                
                if err_lcp <= tol_lcp
                    solved = 1; ite_lcp{i}(j,7) = count;
                end
                if err_lcp < best_err
                    best_err = err_lcp; best_count = count;
                end
                count = count+1;
            end
            ite_lcp{i}(j,8:9) = [best_err, best_count];
        end
        
        fprintf('for the solver %s, the minimal coefficient of the Delassus part is %.4e\n',...
            Algo_name_matlab{i}, mean( abs(min_coef) ) );
    end
    
    fprintf('\n');
    
    for i=1:s_algo_matlab
        SOL_after_perturb = ite_lcp{i}(:,[1 4 7]);
        if size(SOL_after_perturb,1) == 0
            fprintf('the solver: %s has no failed\n', Algo_name_matlab{i} );
            fprintf('\n');
            continue;
        end
        for j=1:3
            nb_success = size(find(SOL_after_perturb(:,j)~=0),1);
            fprintf('using the %dth perturb, the %s successes in %.4g%% of cases\n',...
                j, Algo_name_matlab{i}, nb_success/size(SOL_after_perturb,1)*100);
        end
        fprintf('\n');
    end
    
    fprintf('\n');
    
    for i=1:s_algo_matlab
        SOL_after_perturb = ite_lcp{i}(:,[1 4 7]);
        nb_ite_in_success = ite_lcp{i}(:,[3 6 9]);
        quality_of_sol = ite_lcp{i}(:,[2 5 8]);
        fprintf('Statistics for the solver: %s\n', Algo_name_matlab{i});
        
        for j=1:3
            fprintf('\n');
            fprintf('Case: %dth perturbation:\n',j);
            Idx_success = SOL_after_perturb(:,j)~=0;
            nb_suc = size(find(Idx_success~=0),1);
            nb_fail = size(find(Idx_success==0),1);
            
            if nb_suc==0
                fprintf('no success\n');
            else
                mean_suc_nbite = mean(nb_ite_in_success(Idx_success,j));
                mean_quality_sol = mean(quality_of_sol(Idx_success,j));
                fprintf('for successes: one need %.4g perturb (in average)\n',...
                    mean_suc_nbite);
                fprintf('with a quality of %.2e\n', mean_quality_sol );
            end
            
            if nb_fail==0
                fprintf('no failed\n');
            else
                mean_failed_nbite = mean(nb_ite_in_success(~Idx_success,j));
                mean_quality_unsol = mean(quality_of_sol(~Idx_success,j));
                fprintf('for failures: one need %.4g perturb (in average)\n',...
                    mean_failed_nbite);
                fprintf('with a quality of %.2e\n', mean_quality_unsol );
            end
        end
        
        fprintf('\n\n');
    end
    
%     fprintf('Results of perturbation on Matlab solvers:\n');
%     count_unsol = zeros(size(List_contact_unsolved));
%     for i=1:size(List_contact_unsolved,1)
%         val = Inf*ones(s_algo_matlab,1);
%         idx = zeros(s_algo_matlab,1);
%         nb_perturb = zeros(s_algo_matlab,1);
%         for k=1:s_algo_matlab
%             Idx_solver = Mat_Err(:,end)==k;
%             Idx_unsol = find(Mat_Err(Idx_solver,2)==List_contact_unsolved(i));
%             assert(size(Idx_unsol,1)<2,'At least 2 similar lines!!');
%             if size(Idx_unsol,1)==1
%                 sol_temp = ite_lcp{k}(Idx_unsol,[1 4 7]);
%                 sol_qual = ite_lcp{k}(Idx_unsol,[2 5 8]);
%                 [val(k),idx(k)] = min(sol_qual);
%                 nb_perturb(k) = sol_temp(idx(k));
%             end
%         end
%         [best_err,idx_perturb] = min(val);
%         best_perturb = idx(idx_perturb);
%         nb_per = nb_perturb(idx_perturb);
%         if nb_per > 0
%             fprintf('The unsolved LCP %d has been solved by %s after %d calls of the %dth perturbation with a lcp_error: %.4g\n',...
%                 List_contact_unsolved(i), Algo_name_matlab{idx_perturb}, ...
%                 nb_per, best_perturb, best_err );
%             count_unsol(i) = List_contact_unsolved(i);
%         else
%             fprintf('The unsolved LCP %d has not been solved, the best solution is found by %s after %d calls of the %dth perturbation with a lcp_error: %.4g\n',...
%                 List_contact_unsolved(i), Algo_name_matlab{idx_perturb}, ...
%                 nb_per, best_perturb, best_err );
%         end
%     end
%     
%     fprintf('\n');
%     if size(List_contact_unsolved,2)~=0
%         lcp_unsol_after_perturb = setdiff(List_contact_unsolved,count_unsol);
%         if size(lcp_unsol_after_perturb,2)==0
%             fprintf('There is no LCP remaining unsolved after perturbations with Maltab solvers\n')
%         else
%             fprintf('There exists LCP remaining unsolved after perturbations:\n')
%             for i=1:size(lcp_unsol_after_perturb,2)
%                 fprintf('%d\n',lcp_unsol_after_perturb(i));
%             end
%         end
%     end
%        
%     fprintf('\n\n');
    

    %% new stats:
    Num_err.inside = zeros(size(nb_LCP_solved,1),1); 
    Num_err.final = Num_err.inside;
    Num_err_uns.inside = zeros(size(nb_LCP_unsolved,1),1); 
    Num_err_uns.final = Num_err_uns.inside;
    
    Ext_val.min = zeros(nb_LCP_solved,1); Ext_val.max = Ext_val.min;
    Ext_val.min_q = Ext_val.min; Ext_val.max_q = Ext_val.min;
    Ext_val.min_uns = zeros(nb_LCP_unsolved,1); Ext_val.max_uns = Ext_val.min_uns;
    Ext_val.min_uns_q = Ext_val.min_uns; Ext_val.max_uns_q = Ext_val.min_uns;

    for i=1:nb_LCP_solved
        A = M{i};
        % Delassus part:
        DM = A(1:3*size(A,1)/4,1:3*size(A,1)/4);
        Ext_val.min(i) = min(DM(:)); Ext_val.max(i) = max(DM(:));
        Ext_val.min_q(i) = min( Q{i}(1:3*size(A,1)/4) );
        Ext_val.max_q(i) = max( Q{i}(1:3*size(A,1)/4) );
        
        dim = size(M{i},1);
        d = ones(dim,1);
        [z,w,err_idx,iter,Geo_info,NE] = lcp_lexico_unique_pivot(M{i}, Q{i}, d);
        Num_err.inside(i) = NE.inside;
        Num_err.final(i) = NE.final;
    end
    
    for i=1:nb_LCP_unsolved
        A = M_uns{i};
        % Delassus part:
        DM = A(1:3*size(A,1)/4,1:3*size(A,1)/4);
        Ext_val.min_uns(i) = min(DM(:)); Ext_val.max_uns(i) = max(DM(:));
        Ext_val.min_uns_q(i) = min( Q_uns{i}(1:3*size(A,1)/4) );
        Ext_val.max_uns_q(i) = max( Q_uns{i}(1:3*size(A,1)/4) );
        
        dim = size(M_uns{i},1);
        d = ones(dim,1);
        [z,w,err_idx,iter,Geo_info,NE] = lcp_lexico_unique_pivot(M_uns{i}, Q_uns{i}, d);
        Num_err_uns.inside(i) = NE.inside;
        Num_err_uns.final(i) = NE.final;
    end
    
    k_figs = 1;
    figure(k_figs)
    histogram(Num_err.inside)
    title('numerical error (solved, minimal)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Num_err.final)
    title('numerical error (solved, at the end)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Num_err_uns.inside)
    title('numerical error (unsolved, minimal)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Num_err_uns.final)
    title('numerical error (unsolved, at the end)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.min)
    title('minimal DM values (solved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.max)
    title('maximal DM values (solved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.min_uns)
    title('minimal DM values (unsolved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.max_uns)
    title('maximal DM values (unsolved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.max_uns./Ext_val.min_uns)
    title('ratio DM values (max/min, unsolved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.min_uns_q)
    title('minimal q values (unsolved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.max_uns_q)
    title('maximal q values (unsolved)')
    
    k_figs = k_figs+1;
    figure(k_figs)
    histogram(Ext_val.max_uns_q./Ext_val.min_uns_q)
    title('ratio q values (max/min, unsolved)')
    
    
    fprintf('TODO: statistics on the contact graph...\n');
    
    %% Is M degenerate or/and sufficient?
    % For M is PSD matrix => M is sufficient
    % M on the FC form (Delassus part only), M is PSD.
    % M on the Anistescu-Potra-Stewart (APS) form, M is not PSD, M is
    % degenerate.
    % So: Is Delassus part degenerate? Is M on APS form sufficient?
    nb_LCP_solved = size(M,1); nb_LCP_unsolved = size(M_uns,1);
    nb_tot = nb_LCP_solved+nb_LCP_unsolved;
    DM = cell(nb_tot,1); degenerate = zeros(nb_tot,1);
    
    % Delassus part:
    for i=1:nb_LCP_solved
        nbc = size(M{i},1)/4;
        assert(floor(nbc)==nbc,'BIGRE!!!')
        s_del = 3*nbc;
        ind_rl = [1:nbc nbc+1:2:s_del];
        
        DM{i} = M{i}(ind_rl,ind_rl);
    end
    for i=1:nb_LCP_unsolved
        nbc = size(M_uns{i},1)/4;
        assert(floor(nbc)==nbc,'BIGRE!!!')
        s_del = 3*nbc;
        ind_rl = [1:nbc nbc+1:2:s_del];
        
        DM{nb_LCP_solved+i} = M_uns{i}(ind_rl,ind_rl);
    end
    
    for i=1:nb_tot
        ns = size(DM{i},1);
        for j=1:min(8,ns)
            subM = DM{i}(1:j,1:j);
            if abs(det(subM))<1e-50
                degenerate(i)=1;
                break;
            end
        end
    end
    
    %% Data recovery: rank, cond, size:
    Data_D = zeros(nb_LCP_solved,3);
    Eig_val = [];
    Data_D_uns = zeros(nb_LCP_unsolved,3);
    Eig_val_uns = [];
    
    for i=1:nb_LCP_solved
        A = M{i};
        
        % Delassus part:
        DM = A(1:3*size(A,1)/4,1:3*size(A,1)/4);
        Data_D(j,1) = size(DM,1);
        Data_D(j,2) = cond(DM);
        Data_D(j,3) = rank(DM);
        Eig_val = [Eig_val;eig(DM)];
    end


    for i=1:nb_LCP_unsolved
        A = M_uns{i};
            
        % Delassus part:
        DM = A(1:3*size(A,1)/4,1:3*size(A,1)/4);
        Data_D_uns(j,1) = size(DM,1);
        Data_D_uns(j,2) = cond(DM);
        Data_D_uns(j,3) = rank(DM);
        Eig_val_uns = [Eig_val_uns;eig(DM)];
    end
      
    %% Results:
    Err_uns = zeros(size(M_uns,1),1);
    for i=1:size(M_uns,1)
        Err_uns(i) = SOL_IterLemke_uns{i,2};
    end
    disp(max(Err_uns))
    disp(size(find(Err_uns>1e-8),1))
    
    msg_err = zeros(size(Err_uns));
    for i=1:size(Err_uns,1)
        msg_err(i) = SOL_IterLemke_uns{i,4};
    end
    vv = find(msg_err~=0);
    disp(vv)
    disp(msg_err(vv))
    
    Err = zeros(size(M,1),1);
    for i=1:size(M,1)
        Err(i) = SOL_IterLemke{i,2};
    end
    disp(max(Err))
    disp(size(find(Err>1e-8),1))
    
    msg_err = zeros(size(Err));
    for i=1:size(Err,1)
        msg_err(i) = SOL_IterLemke{i,4};
    end
    vv = find(msg_err~=0);
    disp(vv)
    disp(msg_err(vv))
    
    figure(1)
    bar(1:size(E_d,1),Data_D(:,1))
    hold on
    bar(1:size(E_d,1),Data_D(:,3))
    
    figure(2)
    is_real = zeros(size(Eig_val));
    for i=1:size(Eig_val,1)
        is_real(i) = isreal(Eig_val(i));
    end
    histogram(Eig_val(is_real))
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Mat_Err contains columns corresponding to:
% [Idx_sol | ind | nc | C1: J^TW^+ | C2: EC | \lambda^T J^TW^+ |
% \beta^T D^T W^+ | \alpha^T (\mu\lambda - H^T \beta) | \lambda | 
% \beta | \alpha | J^TW^+ | D^T W^+ | \mu\lambda - H^T \beta | Idx_max,
% Idx_solver],
% with:
% 
% Idx_sol = 0 for unsolved and 1 for solved LCP.
%
% ind: numerous of contact
%
% nc: number of contacts: n_mat
%
% C1: Relative Normal Velocity Criteria (linked with tolerance 1).
%
% C2: Kinetic Energy Criteria (linked with tolerance 2).
% 
% From 6 to 14: number of failure in all data related to the number 
% of contact (4*3 data). The last one corresponds to Impulse outside to the
% the Coulomb's Cone.
%
% Idx_max: mxn matrices: indication on the origin (Idx_mex(1,:)) of the 
% maximal part of the error and it value (Idx_mex(2,:)).
%          mod(Idx_max(1,:),10)= 1: coming from z^T.W
%                                2: coming from z<0
%                                3: coming from W<0
% 
% Idx_solver: number corresponding to the solver: 1 for Lemke, 2 for
% IterLemke and 3 for lexicoLemkeoptplus. (0 for solvers from C++).
%
function [Mat,Idx_max_val] = Calc_Mat_Err(Idx_sol,Mat,Idx_max_val,z,W,...
    nc,M,tol_p,ind,Idx_solver)

    if M~=-1 % no information on mass
        fprintf('Be careful: the error does not take care about the floe mass\n');
    end

    [part,nZ,nW,e] = lcp_error(z,W);
    
    tol = min(tol_p(1),tol_p(2));
    
    tol_I = tol_p(2)/tol_p(1); % for impulse
    
    Mat_temp = [Idx_sol ind nc zeros(1,11) 0 0 0 Idx_solver];
    k = size(Idx_max_val,2); bool_fulfill = 1;
    
    if e >= tol
        % Relative Normal Velocity Criteria:
        W_N = W(1:nc);
        if size(find(W_N<0),1)~=0
            Mat_temp(4) = sum(W_N(W_N<0));
        end
        
        % Kinetic Energy Criteria:
        prod = z.*W;
        prod = prod(1:3*nc);
        if size(find(prod>0),1)~=0
            Mat_temp(5) = sum(prod(prod>0));
        end
        
        % Checking all element of part:
        for i=1:nc
            if part(i)>tol_p(2)
                Mat_temp(6) = Mat_temp(6)+1;
            end
        end
        for i=nc+1:3*nc
            if part(i)>tol_p(2)
                Mat_temp(7) = Mat_temp(7)+1;
            end
        end
        for i=3*nc+1:4*nc
            if part(i)>tol_p(2)
                Mat_temp(8) = Mat_temp(8)+1;
            end
        end
        for i=4*nc+1:4*nc+nZ
            if part(i)>tol_I
                ncz = find(z==-part(i));
                if ncz(1)<=nc
                    Mat_temp(9) = Mat_temp(9)+1;
                elseif ncz(1)<=3*nc
                    Mat_temp(10) = Mat_temp(10)+1;
                else
                    Mat_temp(11) = Mat_temp(11)+1;
                end
            end
        end
        for i=4*nc+nZ+1:4*nc+nZ+nW
            if part(i)>tol_p(1)
                ncw = find(W==-part(i));
                if ncw(1)<=nc
                    Mat_temp(12) = Mat_temp(12)+1;
                elseif ncw(1)<=3*nc
                    Mat_temp(13) = Mat_temp(13)+1;
                else
                    Mat_temp(14) = Mat_temp(14)+1;
                end
            end
        end
        
        % information on the origin of the max part of the error
        err_temp = 0; part_temp = part; Idx_val = [];
        while err_temp<=e/2
            [val,idx] = max(part_temp);
            err_temp = err_temp+val;
            part_temp(idx) = 0;
            
            if idx <= 4*nc
                Idx_max = idx*10+1;
            elseif idx <= 4*nc+nZ
                ncz = find(z==-val);
                Idx_max = ncz(1)*10+2;
            else
                ncw = find(W==-val);
                Idx_max = ncw(1)*10+3;
            end
            Idx_val = [Idx_val [Idx_max;val] ];
        end
        Mat_temp(15:17) = [size(Idx_val,2) Idx_val(1,1) Idx_val(2,1)];
        
    else
        assert(all(part<tol),'BIGRE!!!');
        Mat_temp = []; bool_fulfill = 0;
    end
    
    Mat = [Mat; Mat_temp];
    if bool_fulfill == 1
        Idx_max_val{k+1} = Idx_val;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [part,nZ,nW,e] = lcp_error(z,W)

    part = [abs(z.*W);abs(z(z<0));abs(W(W<0))];
    nZ = size(z(z<0),1);
    nW = size(W(W<0),1);
    e = sum(part);     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = random_perturb(alpha,M)

    min_coef = Inf;

    s_del = 3*size(M,1)/4;
    assert(floor(s_del)==s_del,'BIGRE!!!');
    for i=1:s_del
        for j=1:s_del
            if (M(i,j)~=0) && (abs(M(i,j)) < min_coef)
                min_coef = M(i,j);
            end
        end
    end
    alpha = alpha*min_coef;

    s = size(M,1);
    for i=1:s
        for j=1:s
            if M(i,j)~=0
                M(i,j) = M(i,j) + alpha*rand;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function M = random_perturb2(alpha,M)

    min_coef = Inf;

    s_del = 3*size(M,1)/4;
    assert(floor(s_del)==s_del,'BIGRE!!!');
    
    for i=1:s_del
        for j=1:s_del
            if (M(i,j)~=0) && (abs(M(i,j)) < min_coef)
                min_coef = M(i,j);
            end
        end
    end
    alpha = alpha*min_coef;
    
    for i=1:s_del
        for j=1:s_del
            if M(i,j)~=0
                M(i,j) = M(i,j) + alpha*rand;
            end
        end
    end     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [M,min_coef] = reduction_I(alpha,M)

    min_coef = Inf;

    s_del = 3*size(M,1)/4;
    assert(floor(s_del)==s_del,'BIGRE!!!');
    s_diff = size(M,1)-s_del;
    
    for i=1:s_del
        for j=1:s_del
            if (M(i,j)~=0) && (abs(M(i,j)) < min_coef)
                min_coef = M(i,j);
            end
        end
    end
    alpha = alpha*min_coef;
    
    Id = [eye(s_del) zeros(s_del,s_diff);...
        zeros(s_diff,s_del) zeros(s_diff,s_diff)];
    M = M+alpha*Id;     
end