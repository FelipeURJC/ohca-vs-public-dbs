function figure3
close all

file_prefix = '../data/selected_features_';
file_sufix = {'4_public_db', '8_public_db', '4_ohca_db', '8_ohca_db'};

% Plot Public and OHCA errobars in the same figure.

file_prefix = '../data/selected_features_';
%file_prefix = '';
file_sufix = {'4','8'};

% STEP 2 : plot data
pos = [0 100 700 250];
close all; figure(1);set(gcf, 'Position', pos);
set(gcf,'PaperPositionMode','auto','PaperOrientation','landscape');

h1=subplot('Position', [0.07 0.1 0.4 0.85]);
h2=subplot('Position', [0.575 0.1 0.4 0.85]);

for l = 1: length(file_sufix), % Each different segment length in a separate figure
    ddbb{1} = [file_sufix{l} '_public_db'];
    ddbb{2} = [file_sufix{l} '_ohca_db'];
    if l == 1; axes(h1); else axes(h2); end;
    for db = 1:2, % Both ddbb for the same segment length in the same figure
        file = [file_prefix ddbb{db}];
        load(file);
        
        Nfeat = length(feats_BST.score);
        
        var_name = {'feats_BST','feats_LLR'};
        alg_name = {'BSTsel'};
%         colores = [1, 0, 0; 0.27, 0.58, 0.29];
        colores = [1, 0, 0; 0.32, 0.70, 0.35];
        th_colores = {'k','b'};
        for a = 1, % Two algorithms: BSTsel and L_1LRsel            
            eval(['feats = ' var_name{a} ';']);
            
            ber = squeeze(feats.bootstrap_metrics(:,:,7));
            mber = mean(ber);
            sber = std(ber);
            
            [minerr, idx] = min(mber);
            threshold = minerr + sber(idx);
            opt_n = find(mber<threshold,1,'last');
            
            %figure(l+2*(a-1));
            hold on
            errorbar(Nfeat:-1:1,mber,sber,'.-','Color',colores(db,:),'Markersize',12,'Linewidth',1); 
            line([0 30], threshold*[1 1], 'Color', [colores(db,:) 0.25],'LineWidth',2); %% Threshold with transparency
            plot(Nfeat-idx+1,minerr,'o','Color',colores(db,:),'Markersize',7,'MarkerFaceColor',colores(db,:),'MarkerEdgeColor','k'); 
            plot(Nfeat-opt_n+1,mber(opt_n),'^','Color',colores(db,:),'Markersize',7,'MarkerFaceColor',colores(db,:),'MArkerEdgeColor','k');
            
            % Vertical lines
            line((Nfeat-idx+1)*[1 1], [0 minerr], 'Color', [colores(db,:) 0.25],'LineWidth',1.5); %% Threshold with transparency
            line((Nfeat-opt_n+1)*[1 1], [0 mber(opt_n)], 'Color', [colores(db,:) 0.25],'LineWidth',1.5); %% Threshold with transparency
            text((Nfeat-opt_n+1),-0.5, sprintf('%d', Nfeat-opt_n+1), 'Color', [colores(db,:) 0.25],'Fontsize',11,'HorizontalAlignment', 'center','Fontweight','bold');
            text((Nfeat-idx+1),-0.5, sprintf('%d', Nfeat-idx+1), 'Color', [colores(db,:) 0.25],'Fontsize',11,'HorizontalAlignment', 'center','Fontweight','bold');
            
            xlabel('Number of features','Fontsize',11)
            ylabel('BER (%)','Fontsize',11)
            
            
            %fig_name = ['./' alg_name{a} '-' file_sufix{l} 's'];
            
            l1 = plot(-10,-10,'.-','Color',colores(1,:),'Markersize',12);
            l2 = plot(-10,-10,'.-','Color',colores(2,:),'Markersize',12);
            
            axis([0.7 30.2 0.5 11])
            legend([l1 l2],{'Public','OHCA'},'FontSize',13);
            legend boxoff
            
            %saveas(gca,fig_name);
            %saveas(gca,[fig_name '.eps']);
        end
    end
    set(gca,'YTick',0:2:12, 'YLim', [0 12]);
    if l == 1
        set(gca,'XTick',0:5:25, 'XLim', [0 30]);
    else
        set(gca,'XTick',[0:5:20 30], 'XLim', [0 30]);
    end
end
