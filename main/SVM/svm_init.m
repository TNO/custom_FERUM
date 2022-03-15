function SV_opt = svm_init(gridsel_option,rbf_tab,varargin)

if gridsel_option == 2
   j  = varargin{1};
   SV = svm(kernel('rbf',rbf_tab(j)));
elseif gridsel_option == 1
   SV = svm(kernel('rbf',5));
elseif gridsel_option == 0
   SV = svm(kernel('rbf',rbf_tab));
end

SV.C                 = Inf;
SV.optimizer         = 'libsvm';
SV.nu                = 0;
SV.alpha_cutoff      = -2;

if gridsel_option == 1
   SV_opt = gridsel(param(SV,'rbf',rbf_tab),{'score=cv;score.folds=3;score.leave_one_out=1'})
else
   SV_opt = SV;
end

