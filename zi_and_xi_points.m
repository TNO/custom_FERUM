function [Z1,Z2,X1,X2,WIP,detJ] = zi_and_xi_points(margi,margj,zmax,nIP)

% Computes z1, z2 and x1, x2 values for Gauss integration - Vectorized version

% Integration limits ( should be -infinity, +infinity on both axes, in theory )
zmin = -zmax;

% Determinant of the jacobian of the transformation between [z1max,z1min]x[z2max,z2min] and [-1,1]x[-1,1]
detJ = (zmax-zmin)^2/4;

% Get integration points and weight in [-1,1], nIP is the number of integration pts
[xIP,wIP] = qrule(nIP);

% Transform integration points coordinates from [-1,1] to [zmax,zmin]
z1 = zmin * ones(size(xIP)) + (zmax-zmin) * ( xIP + ones(size(xIP)) ) / 2;
z2 = z1;

% Transform z1 to x1 ( makes use of u_to_x.m function )
i = 1;
z = z1;
marg = margi;
switch marg(i,1)
   case 1  % Normal distribution
      x(i,:) = z(i,:) * marg(i,3) + marg(i,2);
   case 2  % Lognormal distribution
      lambda = marg(i,5);
      zeta = marg(i,6);
      x(i,:) = exp( z(i,:) * zeta + lambda  );
   case 3  % Gamma distribution
      lambda = marg(i,5);
      k = marg(i,6);
      mean = marg(i,2);
      for j = 1 : nIP
         normal_val = normcdf(z(i,j));
         x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
      end
   case 4  % Shifted exponential distribution
      lambda = marg(i,5);
      x_zero = marg(i,6);
      x(i,:) = x_zero + 1/lambda * log( 1 ./ ( 1 - normcdf(z(i,:)) ) );
   case 5  % Shifted Rayleigh distribution
      a = marg(i,5);
      x_zero = marg(i,6);
      x(i,:)= x_zero + a * ( 2*log( 1 ./ (1-normcdf(1,z(i,:))) ) ) .^0.5;
   case 6  % Uniform distribution
      a = marg(i,5);
      b = marg(i,6);
      x(i,:) = a + (b-a) * normcdf(1,z(i,:));
   case 7  % Beta distribution
      q = marg(i,5);
      r = marg(i,6);
      a = marg(i,7);
      b = marg(i,8);
      mean = marg(i,2);
      for j = 1 : nIP
         normal_val = normcdf(z(i,j));
         x01 = fminbnd('zero_beta',0,1,optimset('fminbnd'),q,r,normal_val);
         % Transform x01 from [0,1] to [a,b] interval
         x(i,j) = a + x01 * ( b - a );
      end
   case 8 % Chi-square distribution
      lambda = 0.5;
      nu = marg(i,5);
      k = nu/2 ;
      mean = marg(i,2);
      for j = 1 : nIP
         normal_val = normcdf(z(i,j));
         x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
      end
   case 11 % Type I largest value distribution ( same as Gumbel distribution )
      u_n = marg(i,5);
      a_n = marg(i,6);
      x(i,:) = u_n - (1/a_n) * log( log( 1 ./ normcdf(z(i,:)) ) );
   case 12 % Type I smallest value distribution
      u_1 = marg(i,5);
      a_1 = marg(i,6);
      x(i,:) = u_1 + (1/a_1) * log( log( 1 ./ ( 1 - normcdf(1,z(i,:)) ) ) );
   case 13 % Type II largest value distribution
      u_n = marg(i,5);
      k = marg(i,6);
      x(i,:) = u_n * log( 1 ./ normcdf(z(i,:)) ) .^ (-1/k);
   case 14 % Type III smallest value distribution
      u_1 = marg(i,5);
      k = marg(i,6);
      epsilon = marg(i,7);
      x(i,:) = epsilon + ( u_1 - epsilon ) * log( 1 ./ ( 1 - normcdf(z(i,:)) ) ) .^(1/k);
   case 15 % Gumbel distribution ( same as type I largest value distribution )
      u_n = marg(i,5);
      a_n = marg(i,6);
      x(i,:) = u_n - (1/a_n) * log( log( 1 ./ normcdf(z(i,:)) ) );
   case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
      u_1 = marg(i,5);
      k = marg(i,6);
      x(i,:) = u_1 * log( 1 ./ ( 1 - normcdf(z(i,:)) ) ) .^(1/k);
   case 18 % (Reserved for Laplace distribution)
   case 19 % (Reserved for Pareto distribution)
   case 31 % vector based custom distribution
      ID      = marg(i,5);
      x(i,:)  = custom_invcdf( z(i,:), ID, 'point');
   case 32
        ID      = marg(i,5);
        location= marg(i,6);
        scale   = marg(i,7);

        x(i,:)  = nonparametric_invcdf(z(i,:), ID, location, scale);
   case 51 % Truncated normal marginal distribution
      mean = marg(i,5);
      stdv = marg(i,6);
      xmin = marg(i,7);
      xmax = marg(i,8);
      x(i,:) = mean + stdv * inv_norm_cdf( ...
         normcdf((xmin-mean)/stdv) + ...
         (normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) * normcdf(z(i,:)) ...
         );
   otherwise
end
x1 = x(i,:);

% Transform z2 to x2 ( makes use of u_to_x.m function )
i = 1;
z = z2;
marg = margj;
switch marg(i,1)
   case 1  % Normal distribution
      x(i,:) = z(i,:) * marg(i,3) + marg(i,2);
   case 2  % Lognormal distribution
      lambda = marg(i,5);
      zeta = marg(i,6);
      x(i,:) = exp( z(i,:) * zeta + lambda  );
   case 3  % Gamma distribution
      lambda = marg(i,5);
      k = marg(i,6);
      mean = marg(i,2);
      for j = 1 : nIP
         normal_val = normcdf(z(i,j));
         x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
      end
   case 4  % Shifted exponential distribution
      lambda = marg(i,5);
      x_zero = marg(i,6);
      x(i,:) = x_zero + 1/lambda * log( 1 ./ ( 1 - normcdf(z(i,:)) ) );
   case 5  % Shifted Rayleigh distribution
      a = marg(i,5);
      x_zero = marg(i,6);
      x(i,:)= x_zero + a * ( 2*log( 1 ./ (1-normcdf(z(i,:))) ) ) .^0.5;
   case 6  % Uniform distribution
      a = marg(i,5);
      b = marg(i,6);
      x(i,:) = a + (b-a) * normcdf(z(i,:));
   case 7  % Beta distribution
      q = marg(i,5);
      r = marg(i,6);
      a = marg(i,7);
      b = marg(i,8);
      mean = marg(i,2);
      for j = 1 : nIP
         normal_val = normcdf(z(i,j));
         x01 = fminbnd('zero_beta',0,1,optimset('fminbnd'),q,r,normal_val);
         % Transform x01 from [0,1] to [a,b] interval
         x(i,j) = a + x01 * ( b - a );
      end
   case 8 % Chi-square distribution
      lambda = 0.5;
      nu = marg(i,5);
      k = nu/2 ;
      mean = marg(i,2);
      for j = 1 : nIP
         normal_val = normcdf(z(i,j));
         x(i,j) = fminsearch('zero_gamma',mean,optimset('fminsearch'),k,lambda,normal_val);
      end
   case 11 % Type I largest value distribution ( same as Gumbel distribution )
      u_n = marg(i,5);
      a_n = marg(i,6);
      x(i,:) = u_n - (1/a_n) * log( log( 1 ./ normcdf(z(i,:)) ) );
   case 12 % Type I smallest value distribution
      u_1 = marg(i,5);
      a_1 = marg(i,6);
      x(i,:) = u_1 + (1/a_1) * log( log( 1 ./ ( 1 - normcdf(z(i,:)) ) ) );
   case 13 % Type II largest value distribution
      u_n = marg(i,5);
      k = marg(i,6);
      x(i,:) = u_n * log( 1 ./ normcdf(z(i,:)) ) .^ (-1/k);
   case 14 % Type III smallest value distribution
      u_1 = marg(i,5);
      k = marg(i,6);
      epsilon = marg(i,7);
      x(i,:) = epsilon + ( u_1 - epsilon ) * log( 1 ./ ( 1 - normcdf(z(i,:)) ) ) .^(1/k);
   case 15 % Gumbel distribution ( same as type I largest value distribution )
      u_n = marg(i,5);
      a_n = marg(i,6);
      x(i,:) = u_n - (1/a_n) * log( log( 1 ./ normcdf(z(i,:)) ) );
   case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
      u_1 = marg(i,5);
      k = marg(i,6);
      x(i,:) = u_1 * log( 1 ./ ( 1 - normcdf(z(i,:)) ) ) .^(1/k);
   case 18 % (Reserved for Laplace distribution)
   case 19 % (Reserved for Pareto distribution)
   case 31 % vector based custom distribution
      ID      = marg(i,5);
      x(i,:)  = custom_invcdf( z(i,:), ID, 'point');    
   case 51 % Truncated normal marginal distribution
      mean = marg(i,5);
      stdv = marg(i,6);
      xmin = marg(i,7);
      xmax = marg(i,8);
      x(i,:) = mean + stdv * inv_norm_cdf( ...
         normcdf((xmin-mean)/stdv) + ...
         (normcdf((xmax-mean)/stdv)-normcdf(1,(xmin-mean)/stdv)) * normcdf(z(i,:)) ...
         );
   otherwise
end
x2 = x(i,:);

Z1  = z1' * ones(1,nIP);
Z2  = ones(nIP,1) * z2;
X1  = x1' * ones(1,nIP);
X2  = ones(nIP,1) * x2;
WIP = wIP' * wIP;