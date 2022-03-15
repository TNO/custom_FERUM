function [ G, grad_g ] = gfunbasic(lsf,x,grad_flag,probdata,analysisopt,gfundata)

% The limit-state function is defined by means of an analytical expression or a Matlab m-function,
% using gfundata(lsf).expression.
% This function is called by gfun.m and evaluates gfundata(lsf).expression.

xname = probdata.name;
nx    = size(x,2);
nrv   = size(x,1);


if isfield(gfundata(lsf),'thetag')
   thetagname = gfundata(lsf).thetagname;
   thetag     = gfundata(lsf).thetag;
   nthetag    = size(gfundata(lsf).thetag,2);
   nthetagv   = size(gfundata(lsf).thetag,1);
else
   nthetag    = 0;
end


if isfield(gfundata(lsf),'cg')
   cgname = gfundata(lsf).cgname;
   cg     = gfundata(lsf).cg;
   ncg    = size(cg,2);
   ncgv   = size(cg,1);
else
   ncg    = 0;
end


expression = gfundata(lsf).expression;


if grad_flag == 'yes'
   gradient_expression = gfundata(lsf).dgdq;
end


if nx > 1

   xc = num2cell(x,2);
   eval(sprintf('[ %s ] = deal(xc{:});',[ char(xname)'; blanks(nrv)]));
   
   if isfield(gfundata(lsf),'thetag')
      thetagc = num2cell(thetag');
      eval(sprintf('[ %s ] = deal(thetagc{:});',[ char(thetagname)'; blanks(nthetagv)]));
   end

   if isfield(gfundata(lsf),'cg')
      cgc = num2cell(cg');
      eval(sprintf('[ %s ] = deal(cgc{:});',[ char(cgname)'; blanks(ncgv)]));
   end

   G = eval(expression);
 
   grad_g = 0;
	
elseif nthetag > 1

   G = zeros(1,nthetag);
	
   xc = num2cell(x');
   eval(sprintf('[ %s ] = deal(xc{:});',[ char(xname)'; blanks(nrv)]));

   if isfield(gfundata(lsf),'cg')
      cgc = num2cell(cg');
      eval(sprintf('[ %s ] = deal(cgc{:});',[ char(cgname)'; blanks(ncgv)]));
   end

   for i = 1:nthetag
      thetagic = num2cell(thetag(:,i)');
      eval(sprintf('[ %s ] = deal(thetagic{:});',[ char(thetagname)'; blanks(nthetagv)]));
      G(i) = eval(expression);
   end

   grad_g = 0;
	
else

   xc = num2cell(x');
   eval(sprintf('[ %s ] = deal(xc{:});',[ char(xname)'; blanks(nrv)]));

   if isfield(gfundata(lsf),'thetag')
      thetagc = num2cell(thetag');
      eval(sprintf('[ %s ] = deal(thetagc{:});',[ char(thetagname)'; blanks(nthetagv)]));
   end

   if isfield(gfundata(lsf),'cg')
      cgc = num2cell(cg');
      eval(sprintf('[ %s ] = deal(cgc{:});',[ char(cgname)'; blanks(ncgv)]));
   end

   G = eval(expression);
	
   if grad_flag == 'yes'
      grad_g = zeros(nrv,1);
      for i = 1:nrv
         grad_g(i) = eval(gradient_expression{i});
      end
   else
      grad_g = 0;
   end

end