function figure3
% Plot boxplots for the four databases (Public-4s, Public-8s, OHCA-4s,
% OHCA-8s), 6 performance metrics (Se, Sp, Pp, ACC, BER, AUC) and five
% algorithms (L1-LR, RF, BAG, BST and SVM).

clear all; close all

%metric_names = {'Se', 'Sp','PP','ACC','ERR','FSC','BER','GME','AUC'};
%metrics_to_plot = [1 2 3 4 7 9];  % Select metrics to plot
metric_names = {'Se (%)', 'Sp (%)','PP','ACC','ERR','FSC','BER (%)','GME','AUC'};
metrics_to_plot = [1 2 7];  % Select metrics to plot

ddbb_names = {'Public-4s', 'Public-8s','OHCA-4s','OHCA-8s'};

% 1. Load results in .mat files
load('./data/paired_bootstrap_4_public.mat', 'alg');
results{1} = alg;
load('./data/paired_bootstrap_4_ohca.mat', 'alg');
results{2} = alg;
load('./data/paired_bootstrap_8_public.mat', 'alg');
results{3} = alg;
load('./data/paired_bootstrap_8_ohca.mat', 'alg');
results{4} = alg;

Na = length(alg);               % Number of algorithms
Nm = length(metrics_to_plot);   % Number of metrics to plot
Nb = length(alg(1).output);     % Number of bootstrap resamples
Nr = length(results);           % Number of databases

% 2. Recast from vector of structs to structs of vectors

for r = 1:Nr,
    alg = results{r};
    for m = 1:Nm,
        metric(r,m).vals = zeros(Nb,Na);
        metric(r,m).name = metric_names{metrics_to_plot(m)};
        for a = 1:Na,
            for b = 1:Nb,
                metric(r,m).vals(b,a) = alg(a).output(b).stats(metrics_to_plot(m));
            end
        end
    end
end

metric = metric([1 3 2 4],:);

% 3. Bloxplots

% v_axis = [83 101; ... Vertical axis for Se subplots
%           93 101; ... Vertical axis for Sp subplots
%           88 101; ... Vertical axis for PP subplots
%           90 101; ... Vertical axis for ACC subplots
%           0  11 ; ... Vertical axis for BER subplots
%           98 100.2];% Vertical axis for AUC subplots
v_axis = [75 100.5; ... Vertical axis for Se subplots
          88 100.5; ... Vertical axis for Sp subplots
          0  13];   % Vertical axis for BER subplots];% Vertical axis for AUC subplots



pos = [10 259 800 500];
figure(1);set(gcf, 'Position', pos, 'PaperPosition',pos);

set(0,'DefaultTextInterpreter','tex');
%set(0,'DefaultLegendInterpreter','tex');

nn = 1; % Subplot index

d_x = 0.19;  pi_x = 0.06; 
d_y = 0.275; pi_y = 0.715;

p_x = pi_x;
p_y = pi_y;
for m = 1:Nm,        
    for r = 1:Nr,
        alg = results{r};
        %subplot(Nm,Nr,nn);

        subplot('Position', [p_x, p_y, d_x, d_y]);
        boxplot(metric(r,m).vals,'labels', {alg(:).name},'symbol','ro','outliersize',3);
        set(gca,'XTickLabel',{' '});
        set(gca,'Fontsize',11);
       
        
        if r == 1 || r == 3
            ylabel(metric_names{metrics_to_plot(m)},'Fontsize', 12,'Fontweight', 'bold');
            if m==2
                set(gca,'YTick',88:2:100);
            end
        else 
             set(gca,'YTickLabel',{' '});
        end
        %if m == 1,
        %    title(ddbb_names{r});
        %end
        if m == Nm,
           set (gca,'TicklabelInterpreter','tex')
            set(gca,'XTickLabels',{'L_1-LR', 'RF', 'BAG', 'BST', 'SVM'},'XTickLabelRotation',60);
        end
        ejes = axis;
        ejes(3:4) = v_axis(m,:);
        axis(ejes);
        
        nn = nn+1;
        
        switch r
            case 2
                p_x = p_x + (d_x + 0.12);
            case 4
                p_x = pi_x;
            otherwise
                p_x = p_x + (d_x + 0.025); 
        end
    end
    p_y = p_y - (d_y + 0.025);
end

text(-2.65, -4.5, '(b) OHCA (4-s left, 8-s right)', 'Fontsize', 12,'Fontweight','bold');
text(-16.25, -4.5, '(a) Public (4-s left, 8-s right)', 'Fontsize', 12,'Fontweight','bold');