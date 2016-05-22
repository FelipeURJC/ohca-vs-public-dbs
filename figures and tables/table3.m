function table3

close all; clc

file_prefix = '../../RESULTS/selected_features_';
file_sufix = {'4_public', '8_public', '4_ohca', '8_ohca'};

ddbb_name = {'PUBLIC-4s', 'PUBLIC-8s', 'OHCA-4s', 'OHCA-8s'};
alg_name = {'BSTsel','L1-LRsel'};
for l = 1:length(file_sufix),
    
    file = [file_prefix file_sufix{l}];
    load(file);
    
    Nfeat = length(feats_BST.score);
    
    var_name = {'feats_BST','feats_LLR'};
    
    for a = 1:2,
        
        eval(['feats = ' var_name{a}]);
        
        ber = squeeze(feats.bootstrap_metrics(:,:,7));
        mber = mean(ber);
        sber = std(ber);
        
        [minerr, idx] = min(mber);
        threshold = minerr + sber(idx);
        opt_n = find(mber<threshold,1,'last');
        Nselected = Nfeat-opt_n+1;
        
        % Select the relevant features for each resample
        preselected = feats.bootstrap_sortedIndex(:, opt_n:Nfeat);
        
        % Score = how many times each feature is selected
        mean_score = hist(preselected(:), 1:Nfeat);
        mean_score = mean_score ./ max(mean_score) * 100;
        
        [sortedScore, sortedFeatures] = sort(mean_score,'descend');
        %feat_index = sortedFeatures(1:Nselected);
        
        figure(l);
        subplot(1,2,a)
        barh(Nfeat:-1:Nfeat-Nselected+1,...
            sortedScore(1:Nselected),...
            'FaceColor',[1,0,0.2]);
        hold on;
        barh(1:Nfeat-Nselected,...
            sortedScore(Nfeat:-1:Nselected+1),...
            'FaceColor',[0.8,0.8,0.8]);
        set(gca,'YTick',1:Nfeat);
        set(gca,'ytickLabel',get_feat_names(sortedFeatures(end:-1:1)),...
            'Fontsize',10)
        xlabel('Feature score','Fontsize',12,'Interpreter','tex');
        axis tight
        leg = legend('Selected features','Discarded features',...
            'Location','southeast');
        set(leg,'Fontsize',12,'Interpreter','Latex');
        legend boxoff
        title_msg = sprintf('%s: %s',ddbb_name{l},alg_name{a});
        title(title_msg,'Fontsize',14,...
            'Interpreter','Latex');               
        
    end
    
end

end

function fn = get_feat_names(n)

names = {'tci',    'tcsc',    'exp',    'expmod', ...
    'cm',    'cvbin',    'frqbin', ...
    'abin' ,   'kurt'  ,  'vfleak',    'M', ...
    'A1' ,   'A2'  ,  'A3' ,   'mav', ...
    'psr' ,   'hilb',    'SamEn' ,   'x1', ...
    'x2'  ,  'x3'  ,  'x4'  ,  'x5', ...
    'bCP' ,   'bWT',    'bW'   , 'Li', ...
    'count1'   , 'count2' ,   'count3' };

fn = names(n);
end