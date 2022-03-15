function SV_opt = svr_init(gridsel_option,rbf_tab,C,epsilon,varargin)

if gridsel_option == 2
   j  = varargin{1};
   SV = svr(kernel('rbf',rbf_tab(j)));
elseif gridsel_option == 1
   SV = svr(kernel('rbf',5))
elseif gridsel_option == 0
   SV = svr(kernel('rbf',rbf_tab));
end

SV.C                 = C;
SV.epsilon           = epsilon;
SV.optimizer         = 'libsvm';
SV.nu                = 0;
SV.use_signed_output = 0;

if gridsel_option == 1
   SV_opt = gridsel(param(SV,'rbf',rbf_tab),{'score=cv;score.folds=3;score.leave_one_out=1;loss="quadratic_loss"'});
else
   SV_opt = SV;
end





