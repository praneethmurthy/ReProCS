%%%Code to generate the comparison of RPCA algorithms for simulated data


clear;
clc;
close all



addpath('comparison/');
addpath('comparison/PROPACK/')
addpath('comparison/private1/')
addpath('comparison/MinMaxSelection/')

%% Data Generation
n = 1000;
t_max = 10000;
t_train = 200;
miss_s = 0;
miss_s_pca = 0;
alpha = 200;
f = 50;
MC = 1;

t_calc = [alpha-1 : alpha : t_max - t_train];

sigma = 1e-5;
theta_degree = 15;

t_reprocs = 0;
t_reprocs_pca = 0;
t_reprocs_pca_off = 0;
t_reprocs_off = 0;
t_stoc_rpca = 0;
t_grasta = 0;
t_ncrpca = 0;
t_gd = 0;


temp_err_L = zeros(MC, t_max - t_train);
temp_err_L_pca = zeros(MC, t_max - t_train);
temp_err_L_pca_off = zeros(MC, t_max - t_train);
temp_err_L_off = zeros(MC, t_max - t_train);
temp_err_L_stoc = zeros(MC, t_max- t_train);
temp_err_L_grasta = zeros(MC, t_max- t_train);
temp_err_L_ncrpca = zeros(MC, t_max- t_train);
temp_err_L_gd = zeros(MC, t_max- t_train);

temp_SE_reprocs = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_reprocs_pca = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_reprocs_pca_off = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_reprocs_off = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_grasta = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_stoc = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_ncrpca = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
temp_SE_gd = zeros(MC, ceil((t_max- t_train)/alpha) - 1);
tic

for mc = 1 : MC
    
    fprintf('Monte-Carlo iteration %d in progress\n', mc);
    
    %%%Generating support set and sparse vectors
    S = zeros(n, t_max);
    T = zeros(n, t_max);
    
    b0 = 0.3;
    rho = 1;
    alpha1 = 100;
    s = 50;
    s_train = s/2;
    alpha1_train = alpha1/2;
    beta = ceil(b0 * alpha1);
    x_max = 10;
    x_min = 5;
%     alpha_train = alpha1;
%     num_changes = floor((t_max - t_train )/beta);
%     
%     num_changes1 = min(floor(alpha1 / beta), ceil(n/s));
%     
%     %training small outlier fraction
% %         beta_train = ceil(.01 * beta);
%     beta_train = 1;
%     num_train = floor(t_train/ beta_train);
%     num_train_1 = min(floor(alpha1_train / beta_train), ceil(n/s_train));
%     fval1 = 0;
%     flag = 0;
%     ii1 = 1;
%     for ii = 1 : num_train
%         if(~flag) %%downward motion
%             if(ii1 <= num_train_1)
%                 bind = fval1 + (ii1 - 1) * s_train/rho + 1;
%                 sind = min(bind - 1 + s_train, n);
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_train_1 + 1)
%                     flag = 1;
%                     ii1 = 1;
%                     fval2 = bind;
%                 end
%             end
%         else
%             if(ii1 <= num_train_1)
%                 bind = max(fval2 - (ii1 - 1) * s_train/rho, 1);
%                 sind = bind - 1 + s_train;
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_train_1 + 1)
%                     flag = 0;
%                     ii1 = 1;
%                 end
%             end
%         end
%         idx = bind : sind;
%         jdx = (ii-1) * beta_train + 1 : ii * beta_train;
%         S(idx, jdx) = x_min + ...
%             (x_max - x_min) * rand(length(idx), beta_train);
%         T(idx, jdx) = 1;
%     end
%     
%     flag = 0;
%     ii1 = 1;
%     fval1 = 0;
%     
%     for ii = 1 : num_changes
%         if(~flag)   %%downward motion
%             if(ii1 <= num_changes1)
%                 bind = fval1 + (ii1 - 1) * s/rho + 1;
%                 sind = min(bind - 1 + s, n);
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_changes1 + 1)
%                     flag = 1;
%                     ii1 = 1;
%                     fval2 = bind;
%                 end
%             end
%         else
%             if(ii1 <= num_changes1)
%                 bind = max(fval2 - (ii1 - 1) * s/rho , 1);
%                 sind = bind - 1 + s;
%                 ii1 = ii1 + 1;
%                 if(ii1 == num_changes1 + 1)
%                     flag = 0;
%                     ii1 = 1;
%                 end
%             end
%         end
%         idx = bind : sind;
%         jdx = t_train + (ii-1) * beta + 1 : t_train + ii * beta;
%         S(idx, jdx) = x_min + ...
%             (x_max - x_min) * rand(length(idx), beta);
%         T(idx, jdx) = 1;
%     end

    rho_train = 0.01;
    rho = 0.1;
    x_min = 10;
    x_max = 20;
    
    BernMat = rand(n, t_max);
    T(:, 1 : t_train) = 1.* (BernMat(:, 1 : t_train) <= rho_train);
    T(:, t_train + 1 : end) = 1 .* (BernMat(:, t_train + 1 : t_max) <= rho);
    S = (x_min + (x_max - x_min) * rand(n, t_max)) .* T;

    
    
    %%%Generate low-rank matrix
    r_0 = 30;
    r_1 = 0;
    r_2 = 0;
    r = r_0 + r_1 + r_2;
    L = zeros(n, t_max);
%     diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0 - r_1), ...
%         ones(1 , r_1)];
    
    diag_entries = [linspace(sqrt(f), sqrt(f)/2, r_0)];
    t_1 = 3000;
    t_2 = 6000;
    P = orth(randn(n, r_0));
    coeff_train = zeros(r_0, t_max);
    
    for cc = 1 : r_0
        coeff_train(cc, :) = -diag_entries(cc) + ...
            2 * diag_entries(cc) * rand(1, t_max);
    end
    
%     theta1 = theta_degree * pi / 180;
%     theta2 = 1.01 * theta1;
%     temp_rot1 = [cos(theta1), sin(theta1); ...
%         -sin(theta1) cos(theta1)];
%     rot_matrix1 = blkdiag(eye(r_0 - 1), temp_rot1);
%     
%     temp_rot2 = [cos(theta2), sin(theta2); ...
%         -sin(theta2) cos(theta2)];
%     rot_matrix2 = blkdiag(eye(r_0 - 1), temp_rot2);
%     
%     Irr = eye(r_0 + r_1);
%     PP1 = P(:, 1 : r_0 + r_1) * rot_matrix1 * Irr(:, 1 : r_0);
%     PP2 = [PP1, P(:, end)] * rot_matrix2 * Irr(:, 1 : r_0);

    B1 = skewdec(n, .5);
    B1 = B1 / norm(B1);
    B2 = skewdec(n, .4);
    B2 = B2 / norm(B2);

    delta1 = 1e-2;
    delta2 = 1 * delta1;
    PP1 = orth(exp(delta1 * B1) * P);
    PP2 = orth(exp(delta2 * B2) * PP1);
    
    L(:, 1:t_1) = P(:, 1:r_0) * coeff_train(:, 1:t_1);
    L(:, t_1+1:t_2) = PP1 * coeff_train(:, t_1+1:t_2);
    L(:, t_2 + 1 : end) = PP2 * coeff_train(:, t_2+1:end);
    M = L + S;
    
    %% Main online robust PCA algorithm section
    
    %%%Algorithm parameters
    K = 8;
    omega = x_min / 2;
    gamma = sqrt(4 * log(n)/n);
    s = ceil((gamma + rho) * n/1.5);
    
    %%%Call to ReProCS
%     ev_thresh = 7.5961e-04;
    ev_thresh = 5e-2;
    fprintf('ReProCS\t');
    tt2 = tic;
    P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-6, 10));
    [L_hat, P_hat, S_hat, T_hat, t_hat, ...
        P_track_full, P_track_new, t_calc] = ...
        AutoReProCS_del_bern2(M(:, t_train + 1 :end),...
        P_init, ev_thresh, alpha, K, omega, s, 5);
    t_reprocs = t_reprocs + toc(tt2);
    
    %%%Call to ReProCS-PCA
%     ev_thresh = 1e-3;
    fprintf('ReProCS-PCA\t');
    tt9 = tic;
    P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-2, 10));
    [L_hat_pca, P_hat_pca, S_hat_pca, T_hat_pca, t_hat_pca, ...
        P_track_full_pca, P_track_new_pca, t_calc_pca] = ...
        ReProCS_PCA(M(:, t_train + 1 :end),...
        P_init, ev_thresh, alpha, K, omega, s, 5);
    t_reprocs_pca = t_reprocs_pca + toc(tt9);

    fprintf('Offline-ReProCS-PCA\t');
    tt11 = tic;
    P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-2, 10));
    [L_hat_pca_off, P_hat_pca_off, S_hat_pca_off, T_hat_pca_off, t_hat_pca_off, ...
        P_track_full_pca_off, P_track_new_pca_off, t_calc_pca_off] = ...
        Offline_ReProCS_PCA(M(:, t_train + 1 :end),...
        P_init, ev_thresh, alpha, K, omega, s, 5);
    t_reprocs_pca_off = t_reprocs_pca_off + toc(tt11);

    
    %%%Call to offline ReProCS
    ev_thresh = 1e-3;
%     gamma = sqrt(5 * log(n) / n);
    fprintf('Offline ReProCS\t');
    tt7 = tic;
    P_init = orth(ncrpca(M(:, 1 : t_train), r_0, 1e-6, 10));
    [L_hat_off, P_hat_off, S_hat_off, T_hat_off, t_hat_off, ...
        P_track_full_off, P_track_new_off] = ...
        AutoReProCS_offline(M(:, t_train + 1 :end),...
        P_init, ev_thresh, alpha, K, omega, s, 5);
    t_reprocs_off = t_reprocs_off + toc(tt7);
    
    %%Call to ORPCA
    fprintf('ORPCA\t');
    tt1 = tic;
    [U_hat_stoc, R_hat_stoc, S_hat_stoc] = stoc_rpca(...
        M(:, t_train + 1 : end), r);
    t_stoc_rpca = t_stoc_rpca + toc(tt1);
    
    %%%Call to GRASTA
    fprintf('GRASTA\t');
    tt3 = tic;
    [L_hat_grasta, S_hat_grasta, U_hat_grasta] = ...
        run_alg_grasta(M(:, t_train + 1:end), r);
    t_grasta = t_grasta + toc(tt3);
    
    %%%Call to AltProj
    fprintf('AltProj\t');
    L_hat_ncrpca = zeros(n, t_max - t_train);
    S_hat_ncrpca = zeros(n, t_max - t_train);
    tt4 = tic;
    cc = 1;
    for ii = 1 : t_max - t_train
        if(~(mod(ii + 1, alpha)))
            %fprintf('ALtMin %d\n', ii);
            [L_hat_ncrpca(:, 1 : ii - t_train), ...
                S_hat_ncrpca(:, 1 : ii- t_train)] = ...
                ncrpca(M(:, t_train + 1 : ii), r_0 + r_1 + r_2, 1e-6, 100);
            P_track_ncrpca{cc} = orth(L_hat_ncrpca(:, 1 : ii- t_train));
            cc = cc + 1;
        end
    end
    t_ncrpca = t_ncrpca + toc(tt4);
    
    params.step_const = 0.5; % step size parameter for gradient descent
    params.max_iter   = 100;  % max number of iterations
    params.tol        = 1e-6;% stop when ||Y-UV'-S||_F/||Y||_F < tol
    
    % alpha_bnd is some safe upper bound on alpha,
    % that is, the fraction of nonzeros in each row of S (can be tuned)
    %gamma = 1.5;
    alpha_bnd = s/n;
    
    %%%Call to Gradient Descent
    fprintf('Gradient Descent\n');
    L_hat_gd = zeros(n, t_max - t_train);
    tt5 = tic;
    cc = 1;
    for ii = 1 : t_max- t_train
        if(~(mod(ii + 1, alpha)))
            %fprintf('GDMin %d\n', ii);
            [U_gd, V_gd] = rpca_gd(M(:, t_train + 1 : ii), r_0 + r_1 + r_2, ...
                alpha_bnd, params);
            L_hat_gd(:, 1 : ii- t_train) = U_gd * V_gd';
            P_track_gd{cc} = orth(U_gd * V_gd');
            cc = cc + 1;
        end
    end
    t_gd = t_gd + toc(tt5);
    S_hat_gd = M(:, t_train + 1  :end) - L_hat_gd;
    
    %%Compute performance metrics
    L_hat_stoc = cell2mat(U_hat_stoc(end)) * R_hat_stoc;
    temp_err_L(mc, :) = ...
        sqrt(mean((L(:, t_train+1:end) - L_hat).^2, 1)) ...
        ./ sqrt(mean(L(:, t_train+1:end).^2, 1));
    temp_err_L_pca(mc, :) = ...
        sqrt(mean((L(:, t_train+1:end) - L_hat_pca).^2, 1)) ...
        ./ sqrt(mean(L(:, t_train+1:end).^2, 1));
    temp_err_L_pca_off(mc, :) = ...
        sqrt(mean((L(:, t_train+1:end) - L_hat_pca_off).^2, 1)) ...
        ./ sqrt(mean(L(:, t_train+1:end).^2, 1));
    temp_err_L_off(mc, :) = ...
        sqrt(mean((L(:, t_train+1:end) - L_hat_off).^2, 1)) ...
        ./ sqrt(mean(L(:, t_train+1:end).^2, 1));
    temp_err_L_grasta(mc, :) = sqrt(mean((L(:, t_train + 1 :end) - ...
        L_hat_grasta).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
    temp_err_L_stoc(mc, :) = sqrt(mean((L(:, t_train + 1 :end) - ...
        L_hat_stoc).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
    temp_err_L_ncrpca(mc, :) = sqrt(mean((L(:, t_train + 1 :end) - ...
        L_hat_ncrpca).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
    temp_err_L_gd(mc, :) = sqrt(mean((L(:, t_train + 1 :end) - ...
        L_hat_gd).^2, 1)) ./ sqrt(mean(L(:, t_train + 1 :end).^2, 1));
    
%     miss_s = ...
%         miss_s + (length(find(S_hat))- length(find(S)))/numel(S);
    miss_s_pca = ...
        miss_s_pca + (length(find(S_hat_pca))- length(find(S)))/numel(S);
    
    
    
    for jj = 1 : length(t_calc)
        P_hat_stoc = cell2mat(U_hat_stoc(max(t_calc(jj), 1)));
        P_hat_grasta = cell2mat(U_hat_grasta(max(t_calc(jj), 1)));
        P_hat_ncrpca = cell2mat(P_track_ncrpca(jj));
        P_hat_gd = cell2mat(P_track_gd(jj));
        if (t_calc(jj) +t_train < t_1)
            temp_SE_reprocs(mc, jj) = ...
                Calc_SubspaceError(P_init, ...
                P(:, 1:r_0));
            temp_SE_reprocs_pca(mc, jj) = ...
                Calc_SubspaceError(P_init, ...
                P(:, 1:r_0));
            temp_SE_reprocs_pca_off(mc, jj) = ...
                Calc_SubspaceError(P_init, ...
                P(:, 1:r_0));
            
            temp_SE_reprocs_off(mc, jj) = ...
                Calc_SubspaceError(P_init, ...
                P(:, 1:r_0));
            temp_SE_grasta(mc, jj) = ...
                Calc_SubspaceError(P_hat_grasta, ...
                P(:, 1 :r_0));
            temp_SE_stoc(mc, jj) = ...
                Calc_SubspaceError(P_hat_stoc, ...
                P(:, 1 :r_0));
            temp_SE_ncrpca(mc, jj) = ...
                Calc_SubspaceError(P_hat_ncrpca, ...
                P(:, 1 :r_0));
            temp_SE_gd(mc, jj) = ...
                Calc_SubspaceError(P_hat_gd, ...
                P(:, 1 :r_0));
        elseif((t_calc(jj) +t_train >= t_1) && (t_calc(jj) + t_train < t_2))
            temp_SE_reprocs(mc, jj) = ...
                Calc_SubspaceError(P_track_full{jj}, PP1);
            temp_SE_reprocs_pca(mc, jj) = ...
                Calc_SubspaceError(P_track_full_pca{jj}, PP1);
            temp_SE_reprocs_pca_off(mc, jj) = ...
                Calc_SubspaceError(P_track_full_pca_off{jj}, PP1);
            temp_SE_reprocs_off(mc, jj) = ...
                Calc_SubspaceError(P_track_full_off{jj}, PP1);
            temp_SE_grasta(mc, jj) = ...
                Calc_SubspaceError(P_hat_grasta,...
                PP1);
            temp_SE_stoc(mc, jj) = ...
                Calc_SubspaceError(P_hat_stoc, PP1);
            temp_SE_ncrpca(mc, jj) = ...
                Calc_SubspaceError(P_hat_ncrpca, PP1);
            temp_SE_gd(mc, jj) = ...
                Calc_SubspaceError(P_hat_gd, PP1);
        else
            temp_SE_reprocs(mc, jj) = ...
                Calc_SubspaceError(P_track_full{jj}, PP2);
            temp_SE_reprocs_pca(mc, jj) = ...
                Calc_SubspaceError(P_track_full_pca{jj}, PP2);
            temp_SE_reprocs_pca_off(mc, jj) = ...
                Calc_SubspaceError(P_track_full_pca_off{jj}, PP2);
            temp_SE_reprocs_off(mc, jj) = ...
                Calc_SubspaceError(P_track_full_off{jj}, PP2);
            temp_SE_grasta(mc, jj) = ...
                Calc_SubspaceError(P_hat_grasta,...
                PP2);
            temp_SE_stoc(mc, jj) = ...
                Calc_SubspaceError(P_hat_stoc, PP2);
            temp_SE_ncrpca(mc, jj) = ...
                Calc_SubspaceError(P_hat_ncrpca, PP2);
            temp_SE_gd(mc, jj) = ...
                Calc_SubspaceError(P_hat_gd, PP2);
        end
    end
end

toc

err_L = mean(temp_err_L, 1);
err_L_pca = mean(temp_err_L_pca, 1);
err_L_pca_off = mean(temp_err_L_pca_off, 1);
err_L_off = mean(temp_err_L_off, 1);
err_L_stoc = mean(temp_err_L_stoc, 1);
err_L_grasta = mean(temp_err_L_grasta, 1);
err_L_ncrpca = mean(temp_err_L_ncrpca, 1);
err_L_gd = mean(temp_err_L_gd, 1);

err_SE_reprocs = mean(temp_SE_reprocs, 1);
err_SE_reprocs_pca = mean(temp_SE_reprocs_pca, 1);
err_SE_reprocs_pca_off = mean(temp_SE_reprocs_pca_off, 1);
err_SE_reprocs_off = mean(temp_SE_reprocs_off, 1);
err_SE_grasta = mean(temp_SE_grasta, 1);
err_SE_stoc = mean(temp_SE_stoc, 1);
err_SE_ncrpca = mean(temp_SE_ncrpca, 1);
err_SE_gd = mean(temp_SE_gd, 1);

t_reprocs = t_reprocs / (MC * t_max)
t_reprocs_pca = t_reprocs_pca / (MC * t_max)
t_reprocs_pca_off = t_reprocs_pca_off / (MC * t_max)
t_reprocs_off = t_reprocs_off / (MC * t_max)
t_stoc_rpca = t_stoc_rpca / (MC * t_max)
t_grasta = t_grasta / (MC * t_max)
t_ncrpca = t_ncrpca / (MC * t_max)
t_gd = t_gd / (MC * t_max)


% figure;
% subplot(311)
% plot(t_calc, log10(SE_Prot_Erot), 'r');
% subplot(312)
% plot(t_calc, log10(SE_Phat_Prot), 'k');
% subplot(313)
% plot(t_calc, log10(SE_Phat_P), 'b');


% figure;
% subplot(211);
% imagesc(S);
% subplot(212);
% imagesc(S_hat_pca);

% figure
% subplot(321);
% imagesc(S(:, t_train + 1 : end));
% subplot(322);
% imagesc(S_hat(:, t_train + 1 : end));
% subplot(323);
% imagesc(S_hat_grasta)
% subplot(324)
% imagesc(S_hat_stoc)
% subplot(325)
% imagesc(S_hat_ncrpca)
% subplot(326)
% imagesc(S_hat_gd)


% figure
% subplot(211);
% plot(err_L(t_train + 1  : end));
% subplot(212);
% plot(err_L_grasta);

stry ='$$\log_{10}\left(\frac{\|\hat{\ell}_t - \ell_t\|}{\|\ell_t\|}\right)$$';
strx = '$$t$$';
figure
plot(1 : alpha : t_max - t_train -1, log10(err_L_grasta(1 : alpha ...
    : t_max - t_train - 1)), 'r*--', 'LineWidth', 4, 'MarkerSize', 13)
hold
plot(2 : alpha : t_max - t_train - 2, log10(err_L_stoc...
    (1 : alpha : t_max - t_train -1)), 'ks--','LineWidth', 4, 'MarkerSize', 13)
plot(1 : alpha : t_max - t_train -1, log10(err_L...
    (1 : alpha : t_max - t_train -1)), 'bo--', 'LineWidth', 4, 'MarkerSize', 13)
%hold
plot(1 : alpha : t_max - t_train -1, log10(err_L_pca...
    (1 : alpha : t_max - t_train -1)), 'rd:', 'LineWidth', 4, 'MarkerSize', 13)
% hold
plot(1 : alpha : t_max - t_train -1, log10(err_L_pca_off_temp...
    (1 : alpha : t_max - t_train -1)), 'rd--', 'LineWidth', 4, 'MarkerSize', 13)
plot(1 : alpha : t_max - t_train -1, log10(err_L_ncrpca...
    (1 : alpha : t_max - t_train -1)), 'g^--', 'LineWidth', 4, 'MarkerSize', 13)
plot(1 : alpha : t_max - t_train -1, log10(err_L_gd...
    (1 : alpha : t_max - t_train -1)), 'mp--', 'LineWidth', 4, 'MarkerSize', 13)
plot(1 : alpha : t_max - t_train -1, log10(err_L_off...
    (1 : alpha : t_max - t_train -1)), 'bh--', 'LineWidth', 4, 'MarkerSize', 13)
legend({'GRASTA', 'ORPCA', 'ReProCS', 'ReProCS-PCA', 'Offline-ReProCS-PCA', 'AltProj', 'GradDesc', 'Offline ReProCS'}, 'FontSize', 18);
% plot(1 : alpha : t_max - t_train -1, log10(err_L_ncrpca_alpha...
%     (1 : alpha : t_max - t_train -1)), 'g^:', 'LineWidth', 4, 'MarkerSize', 13)
% plot(1 : alpha : t_max - t_train -1, log10(err_L_rpca_alpha...
%     (1 : alpha : t_max - t_train -1)), 'rd:', 'LineWidth', 4, 'MarkerSize', 13)

xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 20);
ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 20);
axis tight
set(gca, 'FontSize', 20)
% legend({'GRASTA', 'ORPCA', 'ReProCS', 'ReProCS-PCA','AltProj', 'GradDesc', 'Offline ReProCS', 'AltProj (\alpha)', 'PCP (\alpha)'}, 'FontSize', 20);

err_L_pca_off_temp = zeros(1, t_max - t_train);
err_L_pca_off_temp(1 : 1600) = err_L_pca_off(1 : 1600);
err_L_pca_off_temp(1601 : 2800) = err_L_pca_off(2801);
err_L_pca_off_temp(2801 : 4599) = err_L_pca_off(2801 : 4599);
err_L_pca_off_temp(5801 : end) = err_L_pca_off(5801 : end);
%err_L_pca_off_temp(5801 : 8699) = err_L_pca_off(5801 : 8699);

% figure;plot(log10(err_L_pca_off_temp))

%  save('data/compare_bern_full_rotation_n1000_mc5_medrop.mat');


stry ='$$\log_{10}\left(\frac{\|\hat{\ell}_t - \ell_t\|}{\|\ell_t\|}\right)$$';
strx = '$$t$$';
figure
plot(t_calc, log10(err_SE_grasta), 'r','LineWidth', 2)
hold
plot(t_calc, log10(err_SE_stoc), 'k','LineWidth', 2)
plot(t_calc, log10(err_SE_reprocs), 'bs', 'LineWidth', 2)
%hold
plot(t_calc, log10(err_SE_reprocs_pca), 'cs', 'LineWidth', 2)
% hold
plot(t_calc, log10(err_SE_reprocs_pca_off_temp), 'rs', 'LineWidth', 2)
plot(t_calc, log10(err_SE_ncrpca), 'g', 'LineWidth', 2)
plot(t_calc, log10(err_SE_gd), 'm', 'LineWidth', 2)
plot(t_calc, log10(err_SE_reprocs_off), 'y', 'LineWidth', 2)
xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 14);
ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 14);
axis tight
xticks
legend('GRASTA', 'ORPCA', 'ReProCS', 'ReProCS-PCA', 'Offline-ReProCS-PCA', 'AltProj', 'GradDesc', 'Offline ReProCS');


% strx = '$$t$$';
% stry ='$$SE(\hat{P}, P)$$';
% figure
% plot(t_calc, log10(err_SE_grasta), 'r','LineWidth', 2)
% hold
% plot(t_calc, log10(err_SE_stoc), 'k','LineWidth', 2)
% plot(t_calc, log10(err_SE_reprocs), 'bs', 'LineWidth', 2)
% % plot(t_calc, log10(err_SE_ncrpca), 'g', 'LineWidth', 2)
% % plot(t_calc, log10(err_SE_gd), 'm', 'LineWidth', 2)
% plot(t_calc, log10(err_SE_reprocs_off), 'y', 'LineWidth', 2)
% xlabel(strx, 'Interpreter', 'LaTeX', 'FontSize', 14);
% ylabel(stry, 'Interpreter', 'LaTeX', 'FontSize', 14);
% axis tight
% legend('GRASTA', 'ORPCA', 'ReProCS', 'AltProj', 'GradDesc', 'Offline ReProCS');