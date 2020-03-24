function [values,variance]=conditional(X,Y)

% returns the conditional values and variance of X given Y and Y given X
% based on Pearson correlation values because if r = 0 data are independent
% (assuming they are jointly normal)
%
% FORMAT [values,variance]=conditional(X,Y)
%
% INPUTS x and y are two vectors of the same length
%
% OUTPUTS values are the conditioned variables X and Y
%         variances are the conditional variances
%
% Cyril Pernet v1 21/05/2012
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012

if size(X)~=size(Y)
    error('X and Y must have the same size')
end

r = corr(X,Y);
Xhat = r*std(X)*Y / std(Y);
Yhat = r*std(Y)*X / std(X);
Cond_stdX = (1-r^2)*std(X);
Cond_stdY = (1-r^2)*std(Y);

values = [Xhat Yhat];
variance = [Cond_stdX^2 Cond_stdY^2];
