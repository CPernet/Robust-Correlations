function [I,value] = madmedianrule(a,type)

% returns a logical vector that flags outliers as 1s based on the MAD estimator
% based on Wilcox 2005 p 101
% 
% FORMAT I = madmedianrule(a,type)
%
% INPUT
% a is a vector or matrix of data
% type indicates the method to use
% for type = 1, MADS = b_n*1.4826*median(abs(a - median(a)) 
% b_n is the finite sample correction factor see William, J Stat Computation
% and Simulation, 81, 11, 2011
% 1.4826 is the consistancy factor (the std) for the Gaussian distribution 
% for type = 2, MADN = median(abs(a - median(a)) ./ 0.6745
% rescaled MAD by the .6745 to estimate the std of the Gaussian
% distribution - See Wilcox 2005 p78
%
% Cyril Pernet / Guillaume Rousselet 
% ---------------------------------
%  Copyright (C) Corr_toolbox 2012

k = 2.2414; % = sqrt(chi2inv(0.975,1)) 
[n,p]=size(a);
M = median(a);
MAD=median(abs(a - repmat(median(a),n,1)));

switch type
    
    case 1
        % Median Absolute Deviation with finite sample correction factor
        if n == 2
            bn=1.197; % 1.196;
        elseif n == 3
            bn=1.49; % 1.495;
        elseif n == 4
            bn=1.36; % 1.363;
        elseif n == 5
            bn=1.217; % 1.206;
        elseif n == 6
            bn=1.189; % 1.200;
        elseif n == 7
            bn=1.138; % 1.140;
        elseif n == 8
            bn=1.127; % 1.129;
        elseif n == 9
            bn=1.101; % 1.107;
        else
            bn=n/(n-0.8);
        end
        
        MADS=repmat((MAD.*1.4826.*bn),n,1);
        I = a > (M+(k.*MADS));
        I = I+isnan(a);
        value = MADS(1,:);
    
    case 2
        % Normalized Median Absolute Deviation
        MADN = repmat((MAD./.6745),n,1); % same as MAD.*1.4826 :-)
        I = (abs(a-repmat(M,n,1)) ./ MADN) > k;
        I = I+isnan(a);
        value = MADN(1,:);
end
