function [nu,x,h,xp,yp]=univar(X)

% Computes the univariate pdf of data and the histogram values.
% Returns the frequency of data per bin (nu), the position of the bins (x) 
% and their size (h). The pdf is returned in yp for the xp values
%
% Cyril Pernet v2 (10-01-2014 - deals with NaN)
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012

for c=1:size(X,2)
    data = X(~isnan(X(:,c)),c);
    mu = mean(data);
    v = var(data);
    
    % get the normal pdf for this distribution
    xp = linspace(min(data),max(data));
    if v <= 0
        error('Variance must be greater than zero')
        return
    end
    arg = ((xp-mu).^2)/(2*v);
    cons = sqrt(2*pi)*sqrt(v);
    yp = (1/cons)*exp(-arg);
    
    % get histogram info using Surges' rule:
    k = round(1 + log2(length(data)));
    [nu{c},x{c}]=hist(data,k);
    h{c} = x{c}(2) - x{c}(1);
end

if c==1; nu = nu{1}; x = x{1}; h = h{1}; end
