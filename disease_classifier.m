close all
clear
clc

%% TRAINING DATA
% read excel sheets from training data
tp_21_table = readtable("training_data/Pooled Data_Oct 5 2023.xlsx", "Sheet", "21 days");
tp_28_table = readtable("training_data/Pooled Data_Oct 5 2023.xlsx", "Sheet", "28 days");
tp_36_table = readtable("training_data/Pooled Data_Oct 5 2023.xlsx", "Sheet", "36 days");
tp_tables = {tp_21_table, tp_28_table, tp_36_table};

% define time points
time_points = ["21 days"; "28 days"; "36 days"];

% values from miRNA ROC analysis
tp_21_ROC_thresh = 0.8224;
tp_28_ROC_thresh = 0.8908;
tp_36_ROC_thresh = 0.7350;
ROC_threshs = [tp_21_ROC_thresh, tp_28_ROC_thresh, tp_36_ROC_thresh];

% data structures
WT_miRNA = {};
KO_miRNA = {};

% for each time point
for i = 1:size(tp_tables,2)
    
    % EXTRACT DATA
    % extract WT data
    WT_miRNA{i} = tp_tables{i}.miRNA(ismember(tp_tables{i}.Genotype, 'WT'),:);
    WT_Myofiber = tp_tables{i}.Myofiber(ismember(tp_tables{i}.Genotype, 'WT'),:);
    WT_index = ismember(tp_tables{i}.Genotype, 'WT');

    % extract KO data
    KO_miRNA{i} = tp_tables{i}.miRNA(ismember(tp_tables{i}.Genotype, 'KO'),:);
    KO_Myofiber = tp_tables{i}.Myofiber(ismember(tp_tables{i}.Genotype, 'KO'),:);  
    KO_index = ismember(tp_tables{i}.Genotype, 'KO');
    
    % CREATE FIGURE
    fig(i) = figure();
    fig(i).Name = time_points(i);

    % MICRO RNA 
    % find threshold # halfway between means 
    % NOT USED
    sep_num_RNA(i) = (mean(WT_miRNA{i}) + mean(KO_miRNA{i}))/2;

    subplot(3,2,[1,3])
    hold on
    plotted = plot(zeros(size(WT_miRNA{i},1)),WT_miRNA{i},'co',...
         zeros(size(KO_miRNA{i},1)),KO_miRNA{i},'rx');
%          [-1, 1], [ROC_threshs(i), ROC_threshs(i)],'k');
    
    % find threshold using LDA
    LDAMdl = fitcdiscr(tp_tables{i}.miRNA, tp_tables{i}.Genotype,'DiscrimType','linear');
    K = LDAMdl.Coeffs(1,2).Const;
    L = LDAMdl.Coeffs(1,2).Linear;
    miRNA_thresh(i) = -K/L(1);
    xlims = xlim;
    h3 = plot(xlims, [miRNA_thresh(i) miRNA_thresh(i)]);
    h3.Color = 'k';

    % add threshold # text
    y_pos = (max([WT_miRNA{i}; KO_miRNA{i}])-...
        min([WT_miRNA{i}; KO_miRNA{i}]))/20 + miRNA_thresh(i);
    text(0.5, y_pos, num2str(miRNA_thresh(i),4))
    
    % axes format
    ylabel('miRNA')
    title('KO vs WT: miRNA') 
    plotted(1).DisplayName = 'WT';
    plotted(end-1).DisplayName = 'KO';
    legend([plotted(1), plotted(end-1)])

    % classifier performance
    fprintf('\n--- %s ---', time_points(i))
    fprintf('\nmiRNA - LDA')
    classifier_performance('miRNA',tp_tables, i, KO_index, WT_index, ...
        miRNA_thresh(i))
%     fprintf('\nmiRNA - ROC')
%     classifier_performance('miRNA',tp_tables, i, KO_index, WT_index, ...
%         ROC_threshs(i))

    % MYOFIBER 
    % find threshold #
    % NOT USED
    sep_num_myo(i) = (mean(WT_Myofiber) + mean(KO_Myofiber))/2;
    
    subplot(3,2,[2,4])
    hold on
    plotted = plot(zeros(size(WT_Myofiber,1)),WT_Myofiber,'co',...
         zeros(size(KO_Myofiber,1)),KO_Myofiber,'rx');
%          [-1, 1], [sep_num_myo(i), sep_num_myo(i)],'k');

    LDAMdl = fitcdiscr(tp_tables{i}.Myofiber, tp_tables{i}.Genotype,'DiscrimType','linear');
    K = LDAMdl.Coeffs(1,2).Const;
    L = LDAMdl.Coeffs(1,2).Linear;
    myo_thresh(i) = -K/L(1);
    xlims = xlim;
    h3 = plot(xlims, [myo_thresh(i) myo_thresh(i)]);
    h3.Color = 'k';

    % add threshold #
    y_pos = (max([WT_Myofiber; KO_Myofiber])-...
        min([WT_Myofiber; KO_Myofiber]))/20 + myo_thresh(i);
    text(0.5, y_pos, num2str(myo_thresh(i),4))
    
    % axes format
    ylabel('myofiber size [um]')
    title('KO vs WT: myofiber size')
    plotted(1).DisplayName = 'WT';
    plotted(end-1).DisplayName = 'KO';
    legend([plotted(1), plotted(end-1)])
    
    % classifier performance
    fprintf('\nmyofiber - LDA')
    classifier_performance('Myofiber',tp_tables, i, KO_index, WT_index, ...
        myo_thresh(i))

    % MYOFIBER vs MICRO RNA
    subplot(3,2,[5, 6])
    plot(WT_Myofiber,WT_miRNA{i},'co',...
         KO_Myofiber,KO_miRNA{i},'rx');
    xlabel('myofiber size [um]')
    ylabel('miRNA [units]')
    title('KO vs WT: miRNA + myofiber size')
    
    % generate linear divider
    LDAMdl = fitcdiscr([tp_tables{i}.Myofiber, tp_tables{i}.miRNA], ...
        tp_tables{i}.Genotype,'DiscrimType','linear');
    draw_boundaries(LDAMdl,'LDA');
    
    legend('WT', 'KO', '')

end

%% TEST BLIND DATA
% read excel sheets from blind data (treated)
% gene therapy
tp_36_GT_table = readtable("blind_data/Pooled Treated Data_Nov 1 2023_ for LDA.xlsx", "Sheet", "Gene Therapy");
% VPA Val proic acid (works through histone modulations)
tp_36_VPA_table = readtable("blind_data/Pooled Treated Data_Nov 1 2023_ for LDA.xlsx", "Sheet", "VPA");

% extract training data
[GT_WT_miRNA, GT_WT_index, GT_WT_T_miRNA, GT_WT_T_index, GT_KO_miRNA, GT_KO_index, GT_KO_T_miRNA, GT_KO_T_index] = read_treated_table(tp_36_GT_table);
[VPA_WT_miRNA, VPA_WT_index, VPA_WT_T_miRNA, VPA_WT_T_index, VPA_KO_miRNA, VPA_KO_index, VPA_KO_T_miRNA, VPA_KO_T_index] = read_treated_table(tp_36_VPA_table);

% add points to figures

% MICRO RNA 
figure()
subplot(1,2,1)
hold on % needs to be after the subplot is called
plotted_WT = plot(zeros(size(WT_miRNA{3},1)),WT_miRNA{3},'co','MarkerSize', 12);
plotted_KO = plot(zeros(size(KO_miRNA{3},1)),KO_miRNA{3},'rx','MarkerSize', 12);
plotted_thresh = plot([-1, 4], [tp_36_ROC_thresh, tp_36_ROC_thresh],'--k');
plotted_GT_WT = plot(ones(size(GT_WT_miRNA,1),1),GT_WT_miRNA, 'co','MarkerSize', 12);
plotted_KO_WT = plot(ones(size(GT_KO_miRNA,1),1),GT_KO_miRNA,'rx','MarkerSize', 12); 
plotted_GT_T_WT = plot(ones(size(GT_WT_T_miRNA,1),1)*2,GT_WT_T_miRNA,'co','MarkerSize', 12); 
plotted_GT_T_KO = plot(ones(size(GT_KO_T_miRNA,1),1)*2,GT_KO_T_miRNA,'rx','MarkerSize', 12); 
plotted_WT(1).DisplayName = 'WT';
plotted_KO(1).DisplayName = 'KO';

% plotted_thresh(1).DisplayName = 'Threshold';
% plotted_GT_WT(1).DisplayName = 'GT WT';
% plotted_GT_T_WT(1).DisplayName = 'GT WT_T';
% plotted_KO_WT(1).DisplayName = 'GT KO';
% plotted_KO_T_WT(1).DisplayName = 'GT KO_T';
% add threshold #
y_pos = (max([WT_miRNA{3}; KO_miRNA{3}; VPA_WT_miRNA; VPA_WT_T_miRNA; VPA_KO_miRNA; VPA_KO_T_miRNA])-...
    min([WT_miRNA{3}; KO_miRNA{3}; VPA_WT_miRNA; VPA_WT_T_miRNA; VPA_KO_miRNA; VPA_KO_T_miRNA]))/20 + sep_num_RNA(i);
text(3, y_pos, num2str(tp_36_ROC_thresh,4))

% axes format
ylabel('miRNA')
title('KO vs WT: miRNA (Gene Therapy Treated)') 
% legend([plotted_WT(1), plotted_KO(1), plotted_thresh(1), plotted_GT_WT(1), plotted_GT_T_WT(1),  plotted_KO_WT(1),  plotted_GT_T_KO(1)])
legend([plotted_WT(1), plotted_KO(1)])
xticks(0:2)
xticklabels({'tp 36', 'GT untreated', 'GT treated'})
xtickangle(45)

subplot(1,2,2)
hold on % needs to be after the subplot is called
plotted_WT = plot(zeros(size(WT_miRNA{3},1)),WT_miRNA{3},'co','MarkerSize', 12);
plotted_KO = plot(zeros(size(KO_miRNA{3},1)),KO_miRNA{3},'rx','MarkerSize', 12);
plotted_thresh = plot([-1, 4], [tp_36_ROC_thresh, tp_36_ROC_thresh],'--k');
plotted_VPA_WT = plot(ones(size(VPA_WT_miRNA,1),1),VPA_WT_miRNA, 'co','MarkerSize', 12);
plotted_KO_WT = plot(ones(size(VPA_KO_miRNA,1),1),VPA_KO_miRNA,'rx','MarkerSize', 12); 
plotted_VPA_T_WT = plot(ones(size(VPA_WT_T_miRNA,1),1)*2,VPA_WT_T_miRNA,'co','MarkerSize', 12); 
plotted_VPA_T_KO = plot(ones(size(VPA_KO_T_miRNA,1),1)*2,VPA_KO_T_miRNA,'rx','MarkerSize', 12); 
plotted_WT(1).DisplayName = 'WT';
plotted_KO(1).DisplayName = 'KO';

% add threshold #
y_pos = (max([WT_miRNA{3}; KO_miRNA{3}; VPA_WT_miRNA; VPA_WT_T_miRNA; VPA_KO_miRNA; VPA_KO_T_miRNA])-...
    min([WT_miRNA{3}; KO_miRNA{3}; VPA_WT_miRNA; VPA_WT_T_miRNA; VPA_KO_miRNA; VPA_KO_T_miRNA]))/20 + sep_num_RNA(i);
text(3, y_pos, num2str(tp_36_ROC_thresh,4))

% axes format
ylabel('miRNA')
title('KO vs WT: miRNA (VPA Treated)') 
% legend([plotted_WT(1), plotted_KO(1), plotted_thresh(1), plotted_GT_WT(1), plotted_GT_T_WT(1),  plotted_KO_WT(1),  plotted_GT_T_KO(1)])
legend([plotted_WT(1), plotted_KO(1)])
xticks(0:2)
xticklabels({'tp 36', 'VPA untreated', 'VPA treated'})
xtickangle(45)

function [WT_miRNA, WT_index, WT_T_miRNA, WT_T_index, KO_miRNA, KO_index, KO_T_miRNA, KO_T_index] = read_treated_table(data_table)
    
    % WT data
    WT_miRNA = data_table.miRNA(ismember(data_table.Genotype, 'WT'),:);
    WT_index = ismember(data_table.Genotype, 'WT');
    
    % treated
    WT_T_miRNA = data_table.miRNA(ismember(data_table.Genotype, 'WT_T'),:);
    WT_T_index = ismember(data_table.Genotype, 'WT_T');
    
    % KO data
    KO_miRNA = data_table.miRNA(ismember(data_table.Genotype, 'KO'),:);
    KO_index = ismember(data_table.Genotype, 'KO');
    
    % treated
    KO_T_miRNA = data_table.miRNA(ismember(data_table.Genotype, 'KO_T'),:);
    KO_T_index = ismember(data_table.Genotype, 'KO_T');

end

