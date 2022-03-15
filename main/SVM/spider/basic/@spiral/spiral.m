
function a = spiral(hyper) 

%=================================================================================
% TOY SPIRAL data generation object
%================================================================================= 
% A=spiral(H) returns a spiral toy data object initialized with hyperparameters H. 
%
% This generates 2 spirals with 4*m points per spiral, where 1*m points are
% points on two perfect spiral and 3*m points are some additive noise on
% the first quarter.
% 
% 
% Hyperparameters, and their defaults
%  n=1           --  a parameter :-) play to explore
%  m=50          --  number of used points per winding
% 
% Model
%
% Methods:
%  generate,train,test
%  Example :  
%   d=gen(spiral({'n=1','m=50'}));
%   [r s0]=train(svm(kernel('rbf',1)),d)
%   plot(s0,[ -20 20 -20 20]);
%=================================================================================
% Reference : 
% Author    : 
% Link      : 
%=================================================================================
  
  a.m=50;
  a.n=2;
  a.noise=1;
  
  
  p=algorithm('spiral');
  a= class(a,'spiral',p);
  
  if nargin==1
    eval_hyper;
  end  
