function  [dmu_dthetaf,dsigma_dthetaf] = dmusigma_dp(marg)

% This function computes the partial derivatives dmu / dpj and dsigma / dpj for j = 1, ..., 4

i = 1;

switch marg(i,1)

   case 1   % Normal distribution

      mean = marg(i,2);
      stdv = marg(i,3);

      dmu_dthetaf.p1 = 1;
      dmu_dthetaf.p2 = 0;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = 0;
      dsigma_dthetaf.p2 = 1;
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 2   % Lognormal distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      lambda = marg(i,5);
      zeta = marg(i,6);

      dmu_dthetaf.p1 = exp(lambda+1/2*zeta^2);
      dmu_dthetaf.p2 = zeta*exp(lambda+1/2*zeta^2);
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = exp(lambda+1/2*zeta^2) * (exp(zeta^2)-1)^(1/2);
      dsigma_dthetaf.p2 = zeta * exp(lambda+1/2*zeta^2) * (exp(zeta^2)-1)^(1/2) + exp(lambda+1/2*zeta^2) / (exp(zeta^2)-1)^(1/2) * zeta * exp(zeta^2);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 3  % Gamma distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      lambda = marg(i,5);
      k = marg(i,6);

      dmu_dthetaf.p1 = -k/lambda^2;
      dmu_dthetaf.p2 = 1/lambda;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = -k^0.5/lambda^2;
      dsigma_dthetaf.p2 = 0.5/(lambda*k^0.5);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 4  % Shifted exponential distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      lambda = marg(i,5);
      x_zero = marg(i,6);

      dmu_dthetaf.p1 = -1/lambda^2;
      dmu_dthetaf.p2 = 1;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = -1/lambda^2;
      dsigma_dthetaf.p2 = 0;
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 5  % Shifted Rayleigh distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      a = marg(i,5);
      x_zero = marg(i,6);

      dmu_dthetaf.p1 = (pi/2)^0.5;
      dmu_dthetaf.p2 = 1;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = (2-pi/2)^0.5;
      dsigma_dthetaf.p2 = 0;
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 6  % Uniform distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      a = marg(i,5);
      b = marg(i,6);

      dmu_dthetaf.p1 = 1/2;
      dmu_dthetaf.p2 = 1/2;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = -1/(2*(3)^0.5);
      dsigma_dthetaf.p2 = 1/(2*(3)^0.5);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 7  % Beta distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      q = marg(i,5);
      r = marg(i,6);
      a = marg(i,7);
      b = marg(i,8);

      dmu_dthetaf.p1 = (b-a)/(q+r)-q*(b-a)/(q+r)^2;
      dmu_dthetaf.p2 = -q*(b-a)/(q+r)^2;
      dmu_dthetaf.p3 = 1-q/(q+r);
      dmu_dthetaf.p4 = q/(q+r);

      dsigma_dthetaf.p1 = -(b-a)/(q+r)^2*(q*r/(q+r+1))^(1/2)+1/2*(b-a)/(q+r)/(q*r/(q+r+1))^(1/2)*(r/(q+r+1)-q*r/(q+r+1)^2);
      dsigma_dthetaf.p2 = -(b-a)/(q+r)^2*(q*r/(q+r+1))^(1/2)+1/2*(b-a)/(q+r)/(q*r/(q+r+1))^(1/2)*(q/(q+r+1)-q*r/(q+r+1)^2);
      dsigma_dthetaf.p3 = -1/(q+r)*(q*r/(q+r+1))^(1/2);
      dsigma_dthetaf.p4 = 1/(q+r)*(q*r/(q+r+1))^(1/2);

   case 8 % Chi-square distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      nu = marg(i,5);
      lambda = 0.5;

      dmu_dthetaf.p1 = 1/2/lambda;
      dmu_dthetaf.p2 = 0;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = 1/4*2^(1/2)/nu^(1/2)/lambda;
      dsigma_dthetaf.p2 = 0;
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 11 % Type I largest value distribution ( same as Gumbel distribution )

      mean = marg(i,2);
      stdv = marg(i,3);
      u_n = marg(i,5);
      a_n = marg(i,6);

      dmu_dthetaf.p1 = 1;
      dmu_dthetaf.p2 = -0.5772156649/a_n^2;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = 0;
      dsigma_dthetaf.p2 = -pi/(a_n^2*6^0.5);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 12 % Type I smallest value distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      u_1 = marg(i,5);
      a_1 = marg(i,6);

      dmu_dthetaf.p1 = 1;
      dmu_dthetaf.p2 = 0.5772156649/a_1^2;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = 0;
      dsigma_dthetaf.p2 = -pi/(a_1^2*6^0.5);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 13 % Type II largest value distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      u_n = marg(i,5);
      k = marg(i,6);

      dmu_dthetaf.p1 = gamma(1-1/k);
      dmu_dthetaf.p2 = u_n*Psi(1-1/k)*gamma(1-1/k)/k^2;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = (gamma(1-2/k)-(gamma(1-1/k))^2)^0.5;
      dsigma_dthetaf.p2 = 1/2*u_n/(gamma(1-2/k)-gamma(1-1/k)^2)^(1/2)*(2*Psi(1-2/k)*gamma(1-2/k)/k^2-2*gamma(1-1/k)^2*Psi(1-1/k)/k^2);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 14 % Type III smallest value distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      u_1 = marg(i,5);
      k = marg(i,6);
      epsilon = marg(i,7);

      dmu_dthetaf.p1 = gamma(1+1/k);
      dmu_dthetaf.p2 = -(u_1-epsilon)*Psi(1+1/k)*gamma(1+1/k)/k^2;
      dmu_dthetaf.p3 = 1-gamma(1+1/k);
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = (gamma(1+2/k)-gamma(1+1/k)^2)^(1/2);
      dsigma_dthetaf.p2 = 1/2*(u_1-epsilon)/(gamma(1+2/k)-gamma(1+1/k)^2)^(1/2)*(-2*Psi(1+2/k)*gamma(1+2/k)/k^2+2*gamma(1+1/k)^2*Psi(1+1/k)/k^2);
      dsigma_dthetaf.p3 = -(gamma(1+2/k)-gamma(1+1/k)^2)^(1/2);
      dsigma_dthetaf.p4 = 0;

   case 15 % Gumbel distribution ( same as type I largest value distribution )

      mean = marg(i,2);
      stdv = marg(i,3);
      u_n = marg(i,5);
      a_n = marg(i,6);

      dmu_dthetaf.p1 = 1;
      dmu_dthetaf.p2 = -0.5772156649/a_n^2;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = 0;
      dsigma_dthetaf.p2 = -pi/(a_n^2*6^0.5);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )

      mean = marg(i,2);
      stdv = marg(i,3);
      u_1 = marg(i,5);
      k = marg(i,6);

      dmu_dthetaf.p1 = gamma(1+1/k);
      dmu_dthetaf.p2 = -u_1*Psi(1+1/k)*gamma(1+1/k)/k^2;
      dmu_dthetaf.p3 = 0;
      dmu_dthetaf.p4 = 0;

      dsigma_dthetaf.p1 = (gamma(1+2/k)-gamma(1+1/k)^2)^(1/2);
      dsigma_dthetaf.p2 = 1/2*u_1/(gamma(1+2/k)-gamma(1+1/k)^2)^(1/2)*(-2*Psi(1+2/k)*gamma(1+2/k)/k^2+2*gamma(1+1/k)^2*Psi(1+1/k)/k^2);
      dsigma_dthetaf.p3 = 0;
      dsigma_dthetaf.p4 = 0;

   case 18 % Laplace distribution

   case 19 % Pareto distribution

   case 51 % Truncated normal marginal distribution

      mean = marg(i,5);
      stdv = marg(i,6);
      xmin = marg(i,7);
      xmax = marg(i,8);

      h = 1/10000;

      dmu_dthetaf.p1 = ( mean_norm_truncated(mean+h,stdv,xmin,xmax) - mean_norm_truncated(mean,stdv,xmin,xmax) ) / h;
      dmu_dthetaf.p2 = ( mean_norm_truncated(mean,stdv+h,xmin,xmax) - mean_norm_truncated(mean,stdv,xmin,xmax) ) / h;
      dmu_dthetaf.p3 = ( mean_norm_truncated(mean,stdv,xmin+h,xmax) - mean_norm_truncated(mean,stdv,xmin,xmax) ) / h;
      dmu_dthetaf.p4 = ( mean_norm_truncated(mean,stdv,xmin,xmax) - mean_norm_truncated(mean,stdv,xmin,xmax-h) ) / h;

      dsigma_dthetaf.p1 = ( stdv_norm_truncated(mean+h,stdv,xmin,xmax) - stdv_norm_truncated(mean,stdv,xmin,xmax) ) / h;
      dsigma_dthetaf.p2 = ( stdv_norm_truncated(mean,stdv+h,xmin,xmax) - stdv_norm_truncated(mean,stdv,xmin,xmax) ) / h;
      dsigma_dthetaf.p3 = ( stdv_norm_truncated(mean,stdv,xmin+h,xmax) - stdv_norm_truncated(mean,stdv,xmin,xmax) ) / h;
      dsigma_dthetaf.p4 = ( stdv_norm_truncated(mean,stdv,xmin,xmax) - stdv_norm_truncated(mean,stdv,xmin,xmax-h) ) / h;


   otherwise

end