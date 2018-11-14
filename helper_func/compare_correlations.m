function stats_results =  compare_correlations(within_ens_corr,between_ens_corr,crit_alpha)
% compare_correlations: computes KS test differences between correlation
% distributions at a given critical alpha value (e.g. 0.05)

for i = 1:length(within_ens_corr)
    true_dat = get_lower_tri(within_ens_corr{i});
    null_dat = get_lower_tri(between_ens_corr{i});
    [stats_results(i).H, stats_results(i).P] = kstest2(true_dat,null_dat,'alpha',crit_alpha,'tail','smaller');
end
    


