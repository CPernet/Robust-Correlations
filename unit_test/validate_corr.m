function validate_corr(name)

% validation routine for input/output format
% also check output values make sense by
% checking against built-in functions if
% available and checking bounds.
%
% Dr Cyril Pernet - University of Edinburgh
% -----------------------------------------
% Copyright (C) Corr_toolbox 2020

%% Pearson
% test p = 1
% check r,t,pval are correct
n=[10 20 40 80 160 320];
for ss = 1:length(n)
    R = mvnrnd([0 0],eye(2),n(ss));
    [r,~,pval,ci] = Pearson(R(:,1),R(:,2),'figure','on');
    [mr,mpval,mcil,mciu] = corrcoef(R(:,1),R(:,2));
    if any(single([r pval ci(1) ci(2)]) ~= single([mr(1,2) mpval(1,2) mcil(1,2) mciu(1,2)]))
        error('r or pvalue don''t match, sample size: %g',ss)
    end
end

