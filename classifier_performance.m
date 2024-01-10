function classifier_performance(measured_quantity,tp_tables, i, KO_index, WT_index, thresh)
    data_table = tp_tables{i}.(measured_quantity);

    % classifier performance
    TP_amount = sum(data_table(KO_index)<thresh);
    P_missed = sum(data_table(KO_index)>thresh);
    TP_rate = TP_amount/sum(KO_index);
    TN_amount = sum(data_table(WT_index)>thresh);
    N_missed = sum(data_table(WT_index)<thresh);
    TN_rate = TN_amount/sum(WT_index);
    
    % print performance
    fprintf('\nthreshold: %.4f', thresh)
    fprintf('\nKOs id''d  : %d', TP_amount)
    fprintf('\nKOs missed: %d', P_missed)
    fprintf('\nWTs id''d  : %d', TN_amount)
    fprintf('\nWTs missed: %d', N_missed)
    fprintf('\n')
end