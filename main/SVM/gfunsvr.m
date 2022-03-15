function [ G, grad_g ] = gfunsvr(lsf,x,grad_flag,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun

% Warning! In FERUM, SVR are trained in the u-space so that calls to gfunsvr.m
% should be done in terms of standard normal vectors, i.e. x = u
% as opposed to any other classical g-functions for which calls are made
% in terms of physical vectors in the x-space.

if isfield(gfundata(lsf),'ng')
   ng = gfundata(lsf).ng;
else
   ng = 1;
end

G = [];
for ig = 1:ng

   gmax            = gfundata(lsf).SVR(ig).gmax;
   gmin            = gfundata(lsf).SVR(ig).gmin;
   hyperparameters = gfundata(lsf).SVR(ig).hyperparameters;
   svm_buf_size    = analysisopt.svm_buf_size;
   gridsel_option  = 0;


   uu = x_to_u(x,probdata);

   G = [ G ; ( gmax - gmin ) * eval_svm(uu,svm_buf_size,hyperparameters,gridsel_option) + gmin ];

end

grad_g = 0;


