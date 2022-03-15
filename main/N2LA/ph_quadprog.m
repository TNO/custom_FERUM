function [H, k, gfundata, rbdo_funvalues_x] = ph_quadprog(thetag,max_iter,rbdo_parameters,rbdo_funvalues_x,rbdo_fundata,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% Polak-He unified method of feasible directions with Armijo Step size
%
% Call Matlab quadprog.m function
%
% solve min max c_k(x) | f_j(x) <= 0,     j = 1 ,...,q, k = 1, ..., p


global rbdo_nform

alpha = rbdo_parameters.alpha;
beta  = rbdo_parameters.beta;
delta = rbdo_parameters.delta;
gamma = rbdo_parameters.gamma;

steplim = rbdo_parameters.steplim;  % Max number of steps in stepsize calculation
iterlim = 50;                       % Max number of iterations in search dir problem for each random outcome of mustart

thetag  = thetag(:);
nthetag = length(thetag);

if isempty(rbdo_parameters.warmstart) == 1
   k = 0;
else
   k = rbdo_parameters.warmstart;
end


p = 1; % p in the ph_quadprog file is equal to one like 'one objective function' returned by fun: c_0 + c_f * p_f
q = rbdo_parameters.problemsize(2); % Number of constraints

if isempty(rbdo_funvalues_x.c) == 1
   [ c, f ] = fun(lsf,thetag,'func',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield);
   phi      = max(c);
   psi      = max(f);
   psi_plus = max(psi,0);
else
   c        = rbdo_funvalues_x.c;
   phi      = rbdo_funvalues_x.phi;
   f        = rbdo_funvalues_x.f;
   psi      = rbdo_funvalues_x.psi;
   psi_plus = max(psi,0);
end

[ dc_dthetag, df_dthetag ] = fun(lsf,thetag,'grad',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield);

A = [ dc_dthetag' df_dthetag' ];
Delta_vec = ones(p,1)*(phi + gamma*psi_plus) - c;
Gamma_vec = ones(q,1)*psi_plus - f;

Hmat = A'*A/delta; % Putting data in quadprog format
f_vec = [ Delta_vec; Gamma_vec ];

LB  = zeros(p+q,1);
UB  = ones(p+q,1);
Aeq = ones(1,p+q);

option = optimset('MaxIter',iterlim,'Display','off','LargeScale','off'); % parameters in quadprog

for i = 1:max_iter,

   itry = 0;
   while 1
      itry = itry + 1;
      mustart = 2*2*(rand(p+q,1)-0.5)/(p+q);  % Empirical setting!
      [ muhat, thetaneg, exitflag ] = quadprog(Hmat,f_vec,[],[],Aeq,1,LB,UB,mustart,option);
      if ~mod(itry,1000), disp(itry), end
      if exitflag==1, break, end
   end

   for kk = 1:p+q  % avoid numerical errors
      if muhat(kk) < 0, muhat(kk) = 0; end
   end

   theta = -0.5*muhat'*Hmat*muhat - f_vec'*muhat;

   H(:,i) = [ thetag' thetag'.*gfundata(lsf).thetag0' ...
      c'      c'.*rbdo_funvalues_x.initial.c0' ...
      phi      max(c'.*rbdo_funvalues_x.initial.c0') ...
      f'      f'.*rbdo_funvalues_x.initial.f0' ...
      psi     max(f'.*rbdo_funvalues_x.initial.f0') ...
      theta ]';

   h = -A*muhat/delta; % Search direction

   % Unified step size calculation (U = unified function)

   [ cnew, fnew ] = fun(lsf, thetag + beta^k*h,'func',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield);
   psinew = max(fnew);
   phinew = max(cnew);
   Unew   = max(phinew - phi - gamma*psi_plus , psinew - psi_plus);

   if Unew > beta^k*alpha*theta

      while Unew > beta^k*alpha*theta

         k = k+1; if (rem(k,10) == 0), fprintf('+'), end
         [ cnew, fnew ] = fun(lsf, thetag + beta^k*h,'func',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield);
         psinew = max(fnew);
         phinew = max(cnew);
         Unew   = max(phinew - phi - gamma*psi_plus , psinew - psi_plus);
         if k > steplim
            disp('WARNING. Upper No. of Steps reaching in Step size calculations')
            break
         end

      end

   else

      while Unew <= beta^k*alpha*theta

         k = k-1; if (rem(k,10) == 0), fprintf('-'), end
         [ cnew, fnew ] = fun(lsf, thetag + beta^k*h,'func',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield);
         psinew = max(fnew);
         phinew = max(cnew);
         Unew   = max(phinew - phi - gamma*psi_plus , psinew - psi_plus);

         if k < -steplim
            disp('WARNING. lower No. of Steps reaching in Step size calculations')
            break
         end

      end

      k = k+1;

   end
   fprintf('\n')

   [ thetag beta^k*h ];
   thetag = thetag + beta^k*h;


   gfundata(lsf).thetag = thetag;

   [ c, f, dc_dthetag, df_dthetag ] = fun(lsf,thetag,'both',rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield); %, ok = 1

   dc_dthetag(1:p,:)   = dc_dthetag(1:p,:) .* ( rbdo_funvalues_x.initial.c0(1:p) * ones(1,nthetag) ) ./ ( ones(p,1) * gfundata(lsf).thetag0' );
   df_dthetag(1:q-1,:) = df_dthetag(1:q-1,:) .* ( rbdo_funvalues_x.initial.f0(1:q-1) * ones(1,nthetag) ) ./ ( ones(q-1,1) * gfundata(lsf).thetag0' );
   df_dthetag(q,:)     = df_dthetag(q,:) * rbdo_funvalues_x.initial.f0(q) ./ gfundata(lsf).thetag0';

   new_c0 = c.*rbdo_funvalues_x.initial.c0;
   new_f0 = abs(f.*rbdo_funvalues_x.initial.f0);

   c = c.*rbdo_funvalues_x.initial.c0./new_c0;
   f = f.*rbdo_funvalues_x.initial.f0./new_f0;

   rbdo_funvalues_x.initial.c0 = new_c0;
   rbdo_funvalues_x.initial.f0 = new_f0;

   phi      = max(c);
   psi      = max(f);
   psi_plus = max(psi,0);

   dc_dthetag(1:p,:)   = dc_dthetag(1:p,:) ./ ( rbdo_funvalues_x.initial.c0(1:p) * ones(1,nthetag) ) .* ( ones(p,1) * gfundata(lsf).thetag0' );
   df_dthetag(1:q-1,:) = df_dthetag(1:q-1,:) ./ ( rbdo_funvalues_x.initial.f0(1:q-1) * ones(1,nthetag) ) .* ( ones(q-1,1) * gfundata(lsf).thetag0' );
   df_dthetag(q,:)     = df_dthetag(q,:) / rbdo_funvalues_x.initial.f0(q) .* gfundata(lsf).thetag0';

   A = [ dc_dthetag' df_dthetag' ];
   Delta_vec = ones(p,1)*(phi + gamma*psi_plus) - c;
   Gamma_vec = ones(q,1)*psi_plus - f;

   Hmat  = A'*A/delta; % Put data in quadprog format
   f_vec = [Delta_vec;Gamma_vec];

   mustart = muhat;

end

H(:,1+max_iter) = [ thetag' thetag'.*gfundata(lsf).thetag0' ...
                    c'      c'.*rbdo_funvalues_x.initial.c0' ...
                    phi      max(c'.*rbdo_funvalues_x.initial.c0') ...
                    f'      f'.*rbdo_funvalues_x.initial.f0' ...
                    psi     max(f'.*rbdo_funvalues_x.initial.f0') ...
                    theta ]';