function [crit_p_value, adj_sig_value] = adjust_p_values(n, alpha, obs_sig)
% put documentaiton
if n <= 40
    [crit_p_value, adj_sig_value] = p_crit_n30(alpha, obs_sig); 

elseif n <= 70
    [crit_p_value, adj_sig_value] = p_crit_n60(alpha, obs_sig); 
    
elseif n <= 100
    [crit_p_value, adj_sig_value] = p_crit_n80(alpha, obs_sig);

elseif n <= 120
    [crit_p_value, adj_sig_value] = p_crit_n100(alpha, obs_sig);
    
else
    crit_p_value = alpha_;
end