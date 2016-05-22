function figure1

fs = 250;

%%%% STEP 2 : plot data
pos = [-1843 259 700 400];
close all; figure(1);set(gcf, 'Position', pos, 'PaperPosition',pos);

%%%% PLOT VF public
load('./data for figures/public_shock'); data=data(1);
s_ecg = data.s_ecg; t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.07 0.76 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', -2:2:2, 'YTickLabel', {'-2.0', '0', '2.0'}, 'YLim', [-2 2], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');

%%%% PLOT NSR public
load('./data for figures/public_noshock'); data=data(1);
s_ecg = data.s_ecg; t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.575 0.76 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', [-0.5 0 1.5], 'YLim', [-0.6 1.5], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');

text(-4.6, -1.5, '(a) ECG from Public DB (VF left, NSR right)', 'Fontsize', 12,'Fontweight','bold');

%%%% PLOT VF ohca
load('./data for figures/ohca_shock'); data=data(3);
s_ecg = data.s_ecg; t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.07 0.375 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', [], 'XLim', [0 8], 'YTick', -2:2:2, 'YLim', [-2 2], 'Fontsize', 11);
ylabel('ECG (mV)');%xlabel('time (s)');
%%%% PLOT VF ohca
load('./data for figures/ohca_shock'); data=data(4);
s_ecg = data.s_ecg; t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.07 0.115 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', -2:2:2, 'YLim', [-2 2], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');

%%%% PLOT PE ohca
load('./data for figures/ohca_noshock'); data=data(1);
s_ecg = 2*data.s_ecg; t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.575 0.375 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', [], 'XLim', [0 8], 'YTick', -0.5:0.5:0.5, 'YLim', [-0.55 0.55], 'Fontsize', 11);
ylabel('ECG (mV)');%xlabel('time (s)');

%%%% PLOT PE ohca
load('./data for figures/ohca_noshock'); data=data(4);
s_ecg = 2*data.s_ecg; t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.575 0.115 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', [-0.5 0 1], 'YLim', [-0.5 1.1], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');


set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition', [1 1 15 15]);


text(-4.3, -1.2, '(b) ECG from OHCA (VF left, PEA right)', 'Fontsize', 12,'Fontweight','bold');