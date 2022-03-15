function [ zero_beta ] = zero_beta(x,q,r,normal_val)

% used for the beta distribution. 
% Notice: definition of gammainc in matlab

%zero_beta = betainc(x,q,r) - normal_val;
zero_beta = abs( betainc(x,q,r) - normal_val );
             