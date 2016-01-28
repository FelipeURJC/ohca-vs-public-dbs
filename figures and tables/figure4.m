function figure4

%%% Load misclassified examples
load('./data for figures/examples_public');
load('./data for figures/examples_OHCA');
%load examples_public.mat
%load examples_OHCA
fs = 250;
t  = (0:size(ECGs,2)-1)/fs;



%%%% STEP 2 : plot data
pos = [-1843 259 700 400];
close all; figure(1);set(gcf, 'Position', pos, 'PaperPosition',pos);

%%%% PLOT VF public
s_ecg=ECGs(1,:);
subplot('Position', [0.07 0.76 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', -1:1:1, 'YTickLabel', {'-1.0', '0', '1.0'}, 'YLim', [-1.2 1.3], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');

%%%% PLOT VF ohca
pos   = strmatch('VF_153', {data.name},'exact'); 
s_ecg = data(pos).ecg;
t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.575 0.76 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', -0.3:0.3:0.3, 'YTickLabel', {'-0.3', '0', '0.3'}, 'YLim', [-0.4 0.4], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');

text(-5.3, -0.75, '(a) Missed VF from public DB (left) and OHCA (right)', 'Fontsize', 12,'Fontweight','bold');

%%%% PLOT NSh public
s_ecg=ECGs(2,:);
t  = (0:size(ECGs,2)-1)/fs;
subplot('Position', [0.07 0.375 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', [], 'XLim', [0 8], 'YTick', [-0.8 0 1.5], 'YTickLabel', {'-0.8', '0', '1.5'}, 'YLim', [-.8 1.6], 'Fontsize', 11);
ylabel('ECG (mV)');%xlabel('time (s)');
%%%% PLOT NSh public
s_ecg=ECGs(3,:);
subplot('Position', [0.07 0.115 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', [-0.6 0 1], 'YTickLabel', {'-0.6', '0', '1.0'}, 'YLim', [-.75 1.2], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');


%%%% PLOT NSh ohca
pos   = strmatch('PE_40', {data.name},'exact'); 
s_ecg = data(pos).ecg;
t = (0:length(s_ecg)-1)/fs;
subplot('Position', [0.575 0.375 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', [], 'XLim', [0 8], 'YTick', [-0.2 0 .2], 'YTickLabel', {'-0.2', '0', '0.2'}, 'YLim', [-.25 .25], 'Fontsize', 11);
ylabel('ECG (mV)');%xlabel('time (s)');

%%%% PLOT NSh ohca
pos   = strmatch('PE_133', {data.name},'exact'); 
s_ecg = data(pos).ecg;
subplot('Position', [0.575 0.115 0.4 0.225]);
plot(t, s_ecg); 
set(gca, 'XTick', 0:2:8, 'XLim', [0 8], 'YTick', [-0.2 0 0.4], 'YTickLabel', {'-0.2', '0', '0.4'}, 'YLim', [-.25 0.5], 'Fontsize', 11);
ylabel('ECG (mV)');xlabel('time (s)');


set(gcf,'PaperPositionMode','manual');
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition', [1 1 15 15]);

text(-6.15, -0.575, '(b) Missed NSh rhythms from public DB (left) and OHCA (right)', 'Fontsize', 12,'Fontweight','bold');