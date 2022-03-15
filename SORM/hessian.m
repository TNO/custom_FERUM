function hess_G = hessian(lsf,dsptx,dsptu,G,probdata,analysisopt,gfundata,femodel,randomfield)


global nfun

transf_type = probdata.transf_type;
marg = probdata.marg;

if isfield(gfundata(lsf),'thetag')
   thetag = gfundata(lsf).thetag;
end
   
switch transf_type
   case 3
      Lo = probdata.Lo;
end

nrv = size(marg,1);
hess_G = zeros(nrv);

% Takes a default value for analysisopt.ffdpara, if not defined in input file
if ~isfield(analysisopt,'ffdpara')
   switch lower(gfundata(lsf).evaluator)
      case 'basic'
         analysisopt.ffdpara = 1000;
      otherwise
         analysisopt.ffdpara = 50;
   end
end
h = 1/analysisopt.ffdpara;

all_perturbed_x_a_step_ahead = zeros(nrv,nrv);
% Size of all_G_a_step_ahead is 1 , nrv

all_perturbed_x_a_step_back = zeros(nrv,nrv);
% Size of all_G_a_step_back is 1 , nrv

all_perturbed_x_both_steps_ahead = zeros(nrv,nrv*(nrv-1)/2);
% Size of all_G_both_steps_ahead is 1 , nrv*(nrv-1)/2)

for i = 1:nrv

   perturbed_u_a_step_ahead = dsptu;
   perturbed_u_a_step_back = dsptu;

   perturbed_u_a_step_ahead(i) = perturbed_u_a_step_ahead(i) + h;
   perturbed_u_a_step_back(i) = perturbed_u_a_step_back(i) - h;
   
   % Transformation from u to x space
   perturbed_x_a_step_ahead = u_to_x(perturbed_u_a_step_ahead,probdata);
   all_perturbed_x_a_step_ahead(:,i) = perturbed_x_a_step_ahead;
   perturbed_x_a_step_back = u_to_x(perturbed_u_a_step_back,probdata);
   all_perturbed_x_a_step_back(:,i) = perturbed_x_a_step_back;

   perturbed_u_a_step_ahead_i = perturbed_u_a_step_ahead;
   
   for j = 1:(i-1)
      
      perturbed_u_both_steps_ahead = perturbed_u_a_step_ahead_i;
      
      perturbed_u_both_steps_ahead(j) = perturbed_u_both_steps_ahead(j) + h;
      
      % Transformation from u to x space
      perturbed_x_both_steps_ahead = u_to_x(perturbed_u_both_steps_ahead,probdata);
      all_perturbed_x_both_steps_ahead(:,(i-2)*(i-1)/2+j) = perturbed_x_both_steps_ahead;
		
   end
      
end

all_x = [ all_perturbed_x_a_step_ahead all_perturbed_x_a_step_back all_perturbed_x_both_steps_ahead ];
[ all_G, dummy ] = gfun(lsf,all_x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);
all_G_a_step_ahead = all_G(1:nrv);
all_G_a_step_back = all_G((nrv+1):(2*nrv));
all_G_both_steps_ahead = all_G((2*nrv+1):(2*nrv+nrv*(nrv-1)/2));

for i = 1:nrv
   hess_G(i,i) = ( all_G_a_step_ahead(i) - 2*G + all_G_a_step_back(i) ) / h^2;
   for j = 1:(i-1)
      hess_G(i,j) = ( ( all_G_both_steps_ahead((i-2)*(i-1)/2+j) - all_G_a_step_ahead(j) ) - ...
                      ( all_G_a_step_ahead(i) - G ) ) / h^2;
      hess_G(j,i) = hess_G(i,j);
   end
end