function varargout = fun(lsf,thetag,grad_flag,rbdo_parameters,rbdo_fundata,rbdo_funvalues_x,probdata,analysisopt,gfundata,femodel,randomfield)


global rbdo_nfun
global rbdo_nform

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   analysisopt.echo_flag = 1;
end
save_echo_flag = analysisopt.echo_flag;
analysisopt.echo_flag = 0;


% Takes a default value for analysisopt.ffdpara_thetag, if not defined in the inputfile
if ~isfield(analysisopt,'ffdpara_thetag')
   switch lower(gfundata(lsf).evaluator)
      case 'basic'
         analysisopt.ffdpara_thetag = 1000;
      otherwise
         analysisopt.ffdpara_thetag = 20;
   end
end

thetagname = gfundata(lsf).thetagname;
nthetag = length(thetag);

cost        = rbdo_fundata.cost;
constraint  = rbdo_fundata.constraint;
p           = rbdo_parameters.problemsize(1);
q           = rbdo_parameters.problemsize(2);
target_beta = rbdo_parameters.target_beta;

if isfield(gfundata(lsf),'cg')
    cgc = num2cell(gfundata(lsf).cg,2);
    eval(sprintf('[ %s ] = deal(cgc{:});',[ char(gfundata(lsf).cgname)'; blanks(length(gfundata(lsf).cg))]));
end

thetagc = num2cell(gfundata(lsf).thetag0 .* thetag,2);
eval(sprintf('[ %s ] = deal(thetagc{:});',[ char(thetagname)'; blanks(nthetag)]));

for i=1:p
   c(i,1) = eval(char(cost(i)));
end

for i=1:(q-1)
   f(i,1) = eval(char(constraint(i)));
end


if strcmp(grad_flag,'grad') | strcmp(grad_flag,'both')

   % Evaluates cost at thetag(i) + allh(i)
   original_thetag = thetag;

   %thetag
   for j=1:nthetag
      if thetag(j) == 0
         allh(j)=1/analysisopt.ffdpara_thetag;
      else
         allh(j) = thetag(j)/analysisopt.ffdpara_thetag;
      end
      thetag(j) = thetag(j) + allh(j);
      thetagc = num2cell(gfundata(lsf).thetag0 .* thetag,2);
      eval(sprintf('[ %s ] = deal(thetagc{:});',[ char(thetagname)'; blanks(nthetag)]));

      for i=1:p
         c_a_step_ahead(i,j) = eval(char(cost(i)));
      end

      for i=1:(q-1)
         f_a_step_ahead(i,j) = eval(char(constraint(i)));
      end
      thetag(j) = original_thetag(j);
   end

   dc_dthetag = (c_a_step_ahead - c * ones(1,nthetag)) ./ ( ones(p,1) * ( gfundata(lsf).thetag0' .* allh ) );
   df_dthetag = (f_a_step_ahead - f * ones(1,nthetag)) ./ ( ones(q-1,1) * ( gfundata(lsf).thetag0' .* allh ) );

end % if strcmp(grad_flag,'grad') | strcmp(grad_flag,'both')


save_flag_sens = gfundata(lsf).flag_sens;
if strcmp(grad_flag,'func')
   gfundata(lsf).flag_sens = 0;
else
   gfundata(lsf).flag_sens = 1;
end


save_thetag = gfundata(lsf).thetag;
gfundata(lsf).thetag = gfundata(lsf).thetag0 .* thetag;
if femodel.save_formresults
   save_cmd = sprintf('save formdata_RBDO_iter=%d.mat probdata analysisopt gfundata femodel;',femodel.RBDO_iter);
   eval(save_cmd);
end
[ formresults, probdata ] = form(lsf,probdata,analysisopt,gfundata,femodel,randomfield);
if femodel.save_formresults
   save_cmd = sprintf('save formresults_RBDO_iter=%d.mat formresults probdata analysisopt gfundata femodel;',femodel.RBDO_iter);
   eval(save_cmd);
end
pf = formresults.pf1;
if strcmp(grad_flag,'grad') | strcmp(grad_flag,'both')
   dbeta_dthetag = formresults.dbeta_dthetag;
end
gfundata(lsf).thetag = save_thetag;
rbdo_nfun = rbdo_nfun + formresults.nfun;
rbdo_nform = rbdo_nform + 1;

gfundata(lsf).flag_sens = save_flag_sens;

beta = - inv_norm_cdf(pf);

if beta == -inf
    beta = -30;
elseif beta == inf
    beta = 30;
end

if strcmp(grad_flag,'func') | strcmp(grad_flag,'both')

   f(q,1) = target_beta - beta;
   
   % Compute the new objective function
   if p == 2
       if isnan(pf)
           disp('WARNING. In fun.m, pf was NaN! pf has been set to zero.')
           pf = 0;
       end
       objective = c(1) + pf * c(2);
   elseif p > 2
       fprintf('ERROR. Too many cost function (i.e. more than 2)!'),pause;
   elseif p < 2
       fprintf('ERROR. Not enough cost function found (i.e. less than 2)!'),pause;
   end
   
   objective = objective ./ rbdo_funvalues_x.initial.c0;

   f(1:q-1,1) = f(1:(q-1),1) ./ rbdo_funvalues_x.initial.f0(1:(q-1),1);

   f(q,1) = f(q,1) / rbdo_funvalues_x.initial.f0(q,1);

end


if strcmp(grad_flag,'grad') | strcmp(grad_flag,'both')

   df_dthetag(q,:) = -dbeta_dthetag';

   % Compute the new objective function gradient
   if p == 2
       if isnan(pf)
           disp('WARNING. In fun.m, pf was NaN! pf has been set to zero.')
           pf = 0;
       end
       dobjective_dthetag = dc_dthetag(1,:) - normpdf(-beta) * dbeta_dthetag(:)' * c(2) + pf * dc_dthetag(2,:);
   elseif p ~= 2
       fprintf('ERROR. Two cost functions expected (i.e. c = c_0 + p_f * c_f)!'),pause;
   end
   
   dobjective_dthetag = dobjective_dthetag ./ ( rbdo_funvalues_x.initial.c0 * ones(1,nthetag) ) .* ( gfundata(lsf).thetag0' );

   df_dthetag(1:q-1,:) = df_dthetag(1:q-1,:) ./ ( rbdo_funvalues_x.initial.f0(1:q-1) * ones(1,nthetag) ) .* ( ones(q-1,1) * gfundata(lsf).thetag0' );

   df_dthetag(q,:) = df_dthetag(q,:) / rbdo_funvalues_x.initial.f0(q) .* gfundata(lsf).thetag0';

end


switch grad_flag
   case 'func'
      varargout{1} = objective;     
      varargout{2} = f;
  case 'grad'
      varargout{1} = dobjective_dthetag;
      varargout{2} = df_dthetag;
  case 'both'
      varargout{1} = objective;
      varargout{2} = f;
      varargout{3} = dobjective_dthetag;
      varargout{4} = df_dthetag;
end


analysisopt.echo_flag = save_echo_flag;