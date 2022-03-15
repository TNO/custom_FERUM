function f = beta_student_cdf(x,Ptarget,nu)
       
% Student's t-distribution (used to find confidence intervals for MC simulations)
% The CDF of the Student's t-distribution with v degrees of freedom is related 
% to the incomplete beta function by:
%       Pr(|X|<x) = betainc(v/(v+x^2),v/2,1/2)
% so
%              |     betainc(v/(v+x^2),v/2,1/2) / 2      for x<0
%       F(x) = |   0.5                                   for x=0
%              | 1 - betainc(v/(v+x^2),v/2,1/2) / 2      for x>0
%

if x > 0
   P = 1 - betainc( nu/(nu+x^2), nu/2, 1/2 ) / 2;
elseif x < 0
   P = betainc( nu/(nu+x^2), nu/2, 1/2 ) / 2;
else
   P = 0.5;
end

f = Ptarget - P;
    
