%%%Demo to implement the automatic ReProCS algorithm for simulated data.
%This folder contains the code accompanying pre-print.
%
%[1] "New Results for Provable Dynamic Robust PCA", Praneeth Narayanamurthy and Namrata Vaswani, arXiv:1705.08948, 2017.
%
%If you use this code please also cite the following papers
%[2] "An online algorithm for separating sparse  and low-dimensional signal sequences from their sum", Han Guo, Chenlu Qiu, and Namrata Vaswani, IEEE Trans. Sig. Proc., 2014.
%[3] "Recursive Robust PCA or Recursive Sparse Recovery in Large but Structure Noise", Chenlu Qiu, Namrata Vaswani, Brain Lois, and Leslie Hogben, IEEE Trans. Info. Theory., 2014.
%[4] "Real-time Robust Principal Components' Pursuit", Chenlu Qiu, and Namrata Vaswani, Allerton, 2010.

%%%If you want to include more than one subspace change, uncomment sections starting with SS2
clear;
clc;
close all


tic
%% Data Generation
n = 1000;
t_max = 5000;
s = ceil(0.05 * n);
t_train = 200;
miss_s = 0;
alpha = 500;
alpha1 = 100;
f = 1;
cnt = 1;
MC = 1;
err_t = zeros(MC, 16);

%% varing initial error and angles
sigmarange = [1e-3, 1e-5, 1e-7, 1e-9];
thetarange = [1, 15, 35, 45];

for ss = 1 : length(sigmarange)
    sigma = sigmarange(ss);
    for ll = 1 : length(thetarange)
        theta_degree = thetarange(ll);
        fprintf('log(Sigma): %d,\t theta: %d\n', ...
            log10(sigma), theta_degree);
        temp_err_L = zeros(MC, t_max - t_train);
        temp_err_SE = zeros(MC, ceil((t_max - t_train)/alpha));

        for mc = 1 : MC
            
            fprintf('Monte-Carlo iteration %d in progress\n', mc);
            
            %%%Generating support set and sparse vectors
            S = zeros(n, t_max);
            rho = 1;
            b0 = 0.1;
            beta = ceil(b0 * alpha1);
            x_max = 25;
            x_min = 10;
            alpha1 = 100;
            num_changes = floor((t_max -t_train)/beta);
            
            num_changes1 = floor(alpha1 / beta);
            
            flag = 0;
            ii1 = 1;
            fval1 = 0;
            for ii = 1 : num_changes
                if(~flag)   %%downward motion
                    if(ii1 <= num_changes1)
                        bind = fval1 + (ii1 - 1) * s/rho + 1;
                        sind = min(bind - 1 + s, n);
                        ii1 = ii1 + 1;
                        if(ii1 == num_changes1 + 1)
                            flag = 1;
                            ii1 = 1;
                            fval2 = bind;
                        end
                    end
                else
                    if(ii1 <= num_changes1)
                        bind = max(fval2 - (ii1 - 1) * s/rho , 1);
                        sind = bind - 1 + s;
                        ii1 = ii1 + 1;
                        if(ii1 == num_changes1 + 1)
                            flag = 0;
                            ii1 = 1;
                        end
                    end
                end
                idx = bind : sind;
                jdx = t_train + (ii-1) * beta + 1 : t_train + ii * beta;
                S(idx, jdx) = x_min + ...
                    (x_max - x_min) * rand(length(idx), beta);
                T(idx, jdx) = 1;
            end
            
            %%%Generate low-rank matrix
            r_0 = 1;
            r_1 = 1;
            r_2 = 1;
		
	    %%SS2	            
            %             r = r_0 + r_1 + r_2;
            r = r_0 + r_1;
            L = zeros(n, t_max);
%             diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0 - r_1), ...
%                 ones(1 , r_1)];
            diag_entries = [sqrt(f) * ones(1, r_0), 1];
            
            
            t_1 = 1000;
	    %%SS2
            %             t_2 = 3000;
            P = orth(randn(n, r));
            coeff_train = zeros(r_0, t_max);
            
            for cc = 1 : r_0
                coeff_train(cc, :) = -diag_entries(cc) + ...
                    2 * diag_entries(cc) * rand(1, t_max);
            end
            
            theta1 = theta_degree * pi / 180;
	    %%SS2
            %             theta2 = 1.01 * theta1;
            temp_rot1 = [cos(theta1), sin(theta1); ...
                -sin(theta1) cos(theta1)];
            rot_matrix1 = blkdiag(eye(r_0 - 1), temp_rot1);
            
	    %%SS2
            %             temp_rot2 = [cos(theta2), sin(theta2); ...
            %                 -sin(theta2) cos(theta2)];
            %             rot_matrix2 = blkdiag(eye(r_0 - 1), temp_rot2);
            
            Irr = eye(r_0 + r_1);
            PP1 = P(:, 1 : r_0 + r_1) * rot_matrix1 * Irr(:, 1 : r_0);
	    %%SS2
            %             PP2 = [PP1, P(:, end)] * rot_matrix2 * Irr(:, 1 : r_0);
            L(:, 1:t_1) = P(:, 1:r_0) * coeff_train(:, 1:t_1);
            %%SS2
	    %             L(:, t_1+1:t_2) = PP1 * coeff_train(:, t_1+1:t_2);
            L(:, t_1 + 1 : end) = PP1 * coeff_train(:, t_1+1:end);
            M = L + S;
            
            %% Main online robust PCA algorithm section
            
            %%%Algorithm parameters
            K = 5;
            omega = x_min/2;
            %jugaad for init errors
            P_init = orth(P(:, 1 : r_0) + sigma * randn(n, r_0)/2.5);
            ev_thresh = 2e-3;
            %             ev_thresh = 0.1/3 * (sin(theta1))^2;
            
            %%%Call to online RPCA function
            [L_hat, P_hat, S_hat, T_hat, t_hat, ...
                P_track_full, P_track_new, t_calc] = ...
                AutoReProCS(M(:, t_train + 1 : end), P_init, ...
                ev_thresh, alpha, K, omega, s, 5);
            
            %%Compute performance metrics
            temp_err_L(mc, :) = ...
                sqrt(mean((L(:, t_train + 1 : end) - L_hat).^2, 1)) ./ ...
                sqrt(mean(L(:, t_train + 1 : end).^2, 1));
            miss_s = ...
                miss_s + (length(find(S_hat))- length(find(S)))/numel(S);
            
            Ea1 = orth(PP1(:, end) - (P_init * (P_init'* PP1(:, end))));
            
            %             err_t(mc, cnt) = t_hat(end) + t_train - t_1;
	    %%SS2
            %             Ea2 = orth(PP2(:, end) - (P_init * (P_init'* PP2(:, end))));

            %%Calculate the subspace error
            for jj = 1 : length(t_calc)
                if (t_calc(jj) + t_train < t_1)
                    temp_SE_Prot_Erot(mc, jj) = ...
                        Calc_SubspaceError([], []);
                    temp_SE_Phat_Prot(mc, jj) = ...
                        Calc_SubspaceError(P_init, []);
                    temp_SE_Phat_P(mc, jj) = ...
                        Calc_SubspaceError(P_init, ...
                        P(:, 1:r_0));

		    %%SS2
                    %                 elseif((t_calc(jj) + t_train >= t_1) && ...
                    %                         (t_calc(jj) +t_train < t_2))
                    %                     temp_SE_Prot_Erot(mc, jj) = ...
                    %                         Calc_SubspaceError(P_track_new{jj}, Ea1);
                    %                     temp_SE_Phat_Prot(mc, jj) = ...
                    %                         Calc_SubspaceError(P_track_full{jj}, PP1(:, end));
                    %                     temp_SE_Phat_P(mc, jj) = ...
                    %                         Calc_SubspaceError(P_track_full{jj}, PP1);

		    %%SS2 replace the t_1 in next line with t_2, and make Ea1 - Ea2 etc
                elseif (t_calc(jj) + t_train >= t_1)
                    temp_SE_Prot_Erot(mc, jj) = ...
                        Calc_SubspaceError(P_track_new{jj}, Ea1);
                    temp_SE_Phat_Prot(mc, jj) = ...
                        Calc_SubspaceError(P_track_full{jj}, PP1(:, end));
                    temp_SE_Phat_P(mc, jj) = ...
                        Calc_SubspaceError(P_track_full{jj}, PP1);
                end
            end
            %fprintf('\n\n');
        end
        err_L(cnt, :) = mean(temp_err_L, 1);
        SE_Prot_Erot(cnt, :) = mean(temp_SE_Prot_Erot, 1);
        SE_Phat_Prot(cnt, :) = mean(temp_SE_Phat_Prot, 1);
        SE_Phat_P(cnt, :) = mean(temp_SE_Phat_P, 1);
        cnt = cnt + 1;
    end
end
toc

%%call to this works only when there are 16 sigma-theta combination, otherwise need to manually generate figures
FigGenReProCS

