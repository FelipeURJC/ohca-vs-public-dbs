function table4

clc;

% Show public 8-s
info = '\t\t\tPUBLIC 8-s\n\t----------------------------------------------\n'; 
filename = '../data/Paired_bootstrap_8_public_db';
metric = 7;
print_table(filename,metric,info)

% Show ohca 8-s
info = '\n\t\t\tOHCA 8-s\n\t-----------------------------------------------\n'; 
filename = '../data/Paired_bootstrap_8_ohca_db';
metric = 7;
print_table(filename,metric,info)


end

function print_table(filename,metric,info)

load(filename) 

metrics_idx_to_plot = metric;
metric_names = {'sen', 'esp','pp','acc','err','fsc','ber','gme'};
alg_names = {'L1-LR', 'RF', 'BAG', 'BST', 'SVM'};

Na = length(alg);               % Number of algorithms
Nm = length(metric_names);      % Number of metrics
[Nb,s] = size(alg(1).output);   % Number of bootstrap resamples and indication of paired bootstrap (s = 2) o tri-paired bootstrap (s = 3)


% 1. Compute paired metrics and recast from vector of structs to structs of vectors
for m = 1:Nm,
    p1metric(m).vals = zeros(Nb,Na); %#ok<*AGROW>
    p1metric(m).name = metric_names{m};
    p2metric(m).vals = zeros(Nb,Na); %#ok<*AGROW>
    p2metric(m).name = metric_names{m};
    p3metric(m).vals = zeros(Nb,Na); %#ok<*AGROW>
    p3metric(m).name = metric_names{m};
    for a = 1:Na,
        for b = 1:Nb,
            p1metric(m).vals(b,a) = alg(a).output(b,1).stats(m) - alg(a).output(b,2).stats(m); % All minus BT
            p2metric(m).vals(b,a) = alg(a).output(b,1).stats(m) - alg(a).output(b,3).stats(m); % All minus LASSO
            p3metric(m).vals(b,a) = alg(a).output(b,2).stats(m) - alg(a).output(b,3).stats(m); % BT minus LASSO
        end
    end
end


% 2. Print table with Decline in BER associated to feature selection for the
%    features selected using the L1-LR or BST algorithms (Table 4)
for m = metrics_idx_to_plot,
    
    level = 95; % In percentaje
    pl    = (100 - level)/2;
    ph    = level - pl;
    
    fprintf(info);
    fprintf('\t ALL-BSTsel \t ALL-L1LRsel \t BSTsel-L1LRsel \n');
    
    for a = 1:Na,
        CIl1 = prctile(p1metric(m).vals(:,a),pl);
        CIh1 = prctile(p1metric(m).vals(:,a),ph);
        CIl2 = prctile(p2metric(m).vals(:,a),pl);
        CIh2 = prctile(p2metric(m).vals(:,a),ph);
        CIl3 = prctile(p3metric(m).vals(:,a),pl);
        CIh3 = prctile(p3metric(m).vals(:,a),ph);
        
        fprintf('%s  \t%2.1f(%2.1f,%2.1f)\t', alg_names{a}, mean(p1metric(m).vals(:,a)), CIl1, CIh1);
        fprintf('%2.1f(%2.1f,%2.1f)\t', mean(p2metric(m).vals(:,a)), CIl2, CIh2);
        fprintf('%2.1f(%2.1f,%2.1f)', mean(p1metric(m).vals(:,a)), CIl3, CIh3);
        
        fprintf('\n');
    end
    
end
end