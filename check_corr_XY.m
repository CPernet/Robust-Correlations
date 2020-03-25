function [X,Y]=check_corr_XY(X,Y)

% routine used by all correlation functions to check the size of X and Y
%
% FORMAT [X,Y]=check_corr_XY(X,Y)
% INPUT   X and Y are the data in of a correlation function
% OUTPUT  X and Y reshaped if necessary, or error message for user
%
% Dr Cyril Pernet - University of Edinburgh
% ----------------------------------------
% Copyright (C) Corr_toolbox 2020


% if X and or Y vectors
if size(X,1) == 1 && size(X,2) > 1; X = X'; end
if size(Y,1) == 1 && size(Y,2) > 1; Y = Y'; end

% if X or Y vector vs. matrix
if size(X,2) == 1 && size(Y,2) > 1; X = repmat(X,1,size(Y,2)); end
if size(Y,1) == 1 && size(X,2) > 1; Y = repmat(Y,1,size(X,2)); end

if any(size(X)~=size(Y)) 
    error('X and Y must have the same size')
end