
function a = adaboost(C,hyper) 
   
%  ADABOOST
%  
%  
% A=ADABOOST(H) implements basic boosting algorithm using max_margin linear hyperplanes. 
%
%
% Hyperparameters, and their defaults
% kmax = 5                -- maximum number of weak learners
%
% Methods:
%  train, test
%
% Model:
%  child                  -- hyperplanes  
% Example:
%
% % Use adaboost with 1-knn as weak learner and validate with 2 fold cross validation.
% c1=[2,0];
% c2=[-2,0];
% X1=randn(50,2)+repmat(c1,50,1);
% X2=randn(50,2)+repmat(c2,50,1);
% 
% [r,a]=train(cv(adaboost(knn),'folds=2'),d);
%
% ========================================================================
% Reference : The boosting approach to machine learning: An overview
% Author    : Robert E. Schapire
% Link      : http://www.cs.princeton.edu/~schapire/uncompress-papers.cgi/msri.ps
% ========================================================================

  
  % model
  if (nargin <1)
      C=dualperceptron('margin=1');
  end
  a.alpha=[];
  a.child={};
  a.kmax=5;
  
  p=algorithm('adaboost');
  a= class(a,'adaboost',p);
 
  if nargin==2,
    eval_hyper;
  end;
  
  a.alpha=ones(a.kmax,1);  
  
  for i=1:a.kmax
  	a.child{i}=C;
  end  
  
  disp([num2str(a.kmax),' (', C.name, ') classifiers '])