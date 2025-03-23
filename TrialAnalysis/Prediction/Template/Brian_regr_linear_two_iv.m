function [test_stat, p_val] = regr_linear_two_iv(dv, iv)

% for stats = regstats(dv, [iv1, iv2]),
% stats.tstat.beta(1) --> intercept
% stats.tstat.beta(2) --> iv1
% stats.tstat.beta(3) --> iv2

iv_pred_error = iv(:,1);
iv_pitch_next = iv(:,2);

stats = regstats(dv, [iv_pred_error, iv_pitch_next], 'linear', 'tstat');

test_stat = stats.tstat.beta(2);
p_val = stats.tstat.pval(2);

end