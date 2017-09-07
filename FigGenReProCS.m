%%Code to generate preliminary figures for ReProCS algo

%clear
%load data/ReProCS_cosamp_varyzz_varytheta_mc30.mat
%addpath('/home/pkurpadn/Downloads/matlab2tikz-master/src/');

str1 = '$$\theta =\ $$';
str3 = '$$\log_{10}SE(\hat{P_0},P_0)=\ $$';
str2 = '$$\log SE(\hat{P}_{\mathrm{rot}}, E_{\mathrm{rot}})$$';
str4 = '$$t$$';

TheInd = [1, 5, 15, 30];
SigInd = [-2, -4, -6, -8];

% t_calc = t_calc(2:end) / 10000;
% SE_Prot_Erot = SE_Prot_Erot(:, 2:end);
% SE_Phat_Prot = SE_Phat_Prot(:, 2:end);
% SE_Phat_P = SE_Phat_P(:, 2:end);


for i = 1 : 4
    figure;
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) + 1, :)), ...
        'r-.', 'LineWidth', 1);
    hold on;
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +2, :)), ...
        'k--', 'LineWidth', 1);
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +3, :)), ...
        'b:', 'LineWidth', 1);
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +4, :)), ...
        'm-', 'LineWidth', 1);
    
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +1, :)), ...
        'rs', 'MarkerFaceColor', [1 0 0], 'LineWidth', 1);    
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +2, :)), ...
        'ks', 'MarkerFaceColor', [0 0 0], 'LineWidth', 1);
    
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +3, :)), ...
        'bs', 'MarkerFaceColor', [0 0 1], 'LineWidth', 1);
    plot(t_calc, log10(SE_Prot_Erot(4* (i - 1) +4, :)), ...
        'ms', 'MarkerFaceColor', 'magenta', 'LineWidth', 1);
    
    h11 = legend('$\theta =\ 1^{\circ}$', '$\theta =\ 5^{\circ}$',...
        '$\theta =\ 15^{\circ}$', '$\theta =\ 30^{\circ}$');
    
%     STRtitle = [str3, num2str(SigInd(i))];
%     title(STRtitle, 'interpreter', 'latex');
%     xlabel(str4, 'interpreter', 'latex');
%     ylabel(str2, 'interpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    set(h11, 'Interpreter', 'latex', 'FontSize', 10);
    %axis tight
    ylim([-10 -1])
    %grid on
%     figname = ['figures/SE_Prot_Erot_zz_', num2str(-SigInd(i)), '.tex'];
%     matlab2tikz(figname, 'width','3.5cm','height','4cm');
%     close;
%     figname = ['figures/SE_Prot_Erot_zz_', num2str(-SigInd(i)), '.eps'];
%     saveas(gca, figname, 'epsc')
%     close
end
%close all;
str2 = '$$\log SE(\hat{P}, P_{\mathrm{rot}})$$';

for i = 1 : 4
    figure;
    plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) + 1, :)), ...
        'r-.', 'LineWidth', 1.5);
    hold on;
    plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +2, :)), ...
        'k--', 'LineWidth', 1.5);
    plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +3, :)), ...
        'b:', 'LineWidth', 1.5);
    plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +4, :)), ...
        'm-', 'LineWidth', 1.5);
    
%     plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +1, :)), ...
%         'rs', 'MarkerFaceColor', [1 0 0], 'LineWidth', 1);
%     plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +2, :)), ...
%         'ks', 'MarkerFaceColor', [0 0 0], 'LineWidth', 1);
%     plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +3, :)), ...
%         'bs', 'MarkerFaceColor', [0 0 1], 'LineWidth', 1);
%     plot(t_calc, log10(SE_Phat_Prot(4* (i - 1) +4, :)), ...
%         'ms', 'MarkerFaceColor', 'magenta', 'LineWidth', 1);
    
    h11 = legend('$\theta =\ 1^{\circ}$', '$\theta =\ 5^{\circ}$',...
        '$\theta =\ 15^{\circ}$', '$\theta =\ 30^{\circ}$');
    
%     STRtitle = [str3, num2str(SigInd(i))];
%     title(STRtitle, 'interpreter', 'latex');
    set(h11, 'Interpreter', 'latex');
%     xlabel(str4, 'interpreter', 'latex');
%     ylabel(str2, 'interpreter', 'latex');
%     yt = get(gca, 'YTick');
%     set(gca, 'FontSize', 16)
%     xt = get(gca, 'XTick');
%     set(gca, 'FontSize', 16)
    axis tight
    ylim([-10 1.1])
%     figname = ['figures/SE_Phat_Prot_zz_', num2str(-SigInd(i)), '.tex'];
%     matlab2tikz(figname, 'width','3.5cm','height','4cm');
%     close;

%     grid on
%     figname = ['figures/SE_Phat_Prot_zz_', num2str(-SigInd(i)), '.eps'];
%     saveas(gca, figname, 'epsc')
%     close
end

str2 = '$$\log SE(\hat{P}, P)$$';
for i = 1 : 4
    figure;
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) + 1, :)), ...
        'r', 'LineWidth', 2);
    hold on;
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +2, :)), ...
        'k', 'LineWidth', 2);
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +3, :)), ...
        'b', 'LineWidth', 2);
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +4, :)), ...
        'g', 'LineWidth', 2);
    
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +1, :)), ...
        'rs', 'MarkerFaceColor', [1 0 0], 'LineWidth', 2);    
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +2, :)), ...
        'ko', 'MarkerFaceColor', [0 0 0], 'LineWidth', 2);
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +3, :)), ...
        'b*', 'MarkerFaceColor', [0 0 1], 'LineWidth', 1);
    plot(t_calc, log10(SE_Phat_P(4* (i - 1) +4, :)), ...
        'g+', 'MarkerFaceColor', [0 1 0], 'LineWidth', 2);
    
    h11 = legend('$\theta =\ 1^{\circ}$', '$\theta =\ 5^{\circ}$',...
        '$\theta =\ 15^{\circ}$', '$\theta =\ 30^{\circ}$');
    
    STRtitle = [str3, num2str(SigInd(i))];
    title(STRtitle, 'interpreter', 'latex');
    set(h11, 'Interpreter', 'latex');
    xlabel(str4, 'interpreter', 'latex');
    ylabel(str2, 'interpreter', 'latex');
    yt = get(gca, 'YTick');
    set(gca, 'FontSize', 16)
    xt = get(gca, 'XTick');
    set(gca, 'FontSize', 16)
    grid on
    axis tight
    
%     axes('position', [.65 .25 .25 .25])
%     box on
%     indexofInterest = (t_calc >= 2000) & (t_calc<= 8000);
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) + 1, (indexofInterest))), ...
%         'r', 'LineWidth', 2);
%     hold on;
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +2, (indexofInterest))), ...
%         'k', 'LineWidth', 2);
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +3, (indexofInterest))), ...
%         'b', 'LineWidth', 2);
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +4, (indexofInterest))), ...
%         'g', 'LineWidth', 2);
%     
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +1, (indexofInterest))), ...
%         'rs', 'MarkerFaceColor', [1 0 0], 'LineWidth', 2);    
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +2, (indexofInterest))), ...
%         'ko', 'MarkerFaceColor', [0 0 0], 'LineWidth', 2);
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +3, (indexofInterest))), ...
%         'b*', 'MarkerFaceColor', [0 0 1], 'LineWidth', 1);
%     plot(t_calc(indexofInterest), ...
%         log10(SE_Phat_P(4* (i - 1) +4, (indexofInterest))), ...
%         'g+', 'MarkerFaceColor', [0 1 0], 'LineWidth', 2);
%     grid on
%     axis tight
%     figname = ['figures/SE_Phat_P_zz_', num2str(-SigInd(i)), '.eps'];
%     saveas(gca, figname, 'epsc')
%     close
end
