function  J_X_Z = dX_dZ(X,Z,marg)

%! NOTE!
%! the same procedure is done in every cases, this whole case stuff could be replaced

i = 1;

switch marg(i,1)
   
   case 1   % Normal distribution
      
      J_Z_X = ones(size(X))/marg(i,3);
      
   case 2   % Lognormal distribution
      
      ksi = sqrt( log( 1 + ( marg(i,3) / marg(i,2) )^2 ) );
      
      J_Z_X = 1 ./ ( ksi * X );
      
   case 3  % Gamma distribution
      
      pdf1 = ferum_pdf(3,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 4  % Shifted exponential distribution
      
      pdf1 = ferum_pdf(4,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 5  % Shifted Rayleigh distribution
      
      pdf1 = ferum_pdf(5,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 6  % Uniform distribution
      
      pdf1 = ferum_pdf(6,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 7  % Beta distribution
      
      pdf1 = ferum_pdf(7,X,marg(i,5),marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 8 % Chi-square distribution
      
      pdf1 = ferum_pdf(8,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X= pdf1./pdf2;
      
   case 11 % Type I largest value distribution ( same as Gumbel distribution )
      
      pdf1 = ferum_pdf(11,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 12 % % Type I smallest value distribution
      
      pdf1 = ferum_pdf(12,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 13 % Type II largest value distribution
      
      pdf1 = ferum_pdf(13,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 14 % Type III smallest value distribution
      
      pdf1 = ferum_pdf(14,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 15 % Gumbel distribution ( same as type I largest value distribution )
      
      pdf1 = ferum_pdf(15,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
      
      pdf1 = ferum_pdf(16,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   case 18 % Laplace distribution
      
   case 19 % Pareto distribution
      
   case 51 % Truncated normal marginal distribution
      
      pdf1 = ferum_pdf(51,X,marg(i,5:8));
      pdf2 = normpdf(Z);
      
      J_Z_X = pdf1./pdf2;
      
   otherwise
      
end

J_X_Z = 1 ./ J_Z_X;