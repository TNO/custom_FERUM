function  J_Z_thetaf = dZ_dthetaf(X,marg)

i = 1;

switch marg(i,1)

   case 1   % Normal distribution
      
      mean = marg(i,2);
      stdv = marg(i,3);

      dzdmean = -1/stdv * ones(size(X));
      dzdstdv = -( X - mean )/(stdv)^2;

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdmean;
      J_Z_thetaf.p2    = dzdstdv;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 2   % Lognormal distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      lambda = marg(i,5);
      zeta = marg(i,6);

      dzdlambda = -1/zeta * ones(size(X));
      dzdzeta = -(log(X)-lambda)/zeta^2;

      dzdmean = dzdlambda * (1/mean + (stdv^2)/(mean^3+(stdv^2)*mean)) + dzdzeta * (-stdv^2)/((mean^3+(stdv^2)*mean)*(log(1+(stdv/mean)^2))^0.5);
      dzdstdv = dzdlambda * (-stdv/(stdv^2+mean^2)) + dzdzeta * stdv/((mean^2+stdv^2)*(log(1+(stdv/mean)^2))^0.5);

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdlambda;
      J_Z_thetaf.p2    = dzdzeta;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 3  % Gamma distribution

      mean   = marg(i,2);
      stdv   = marg(i,3);
      lambda = marg(i,5);
      k      = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(3,X,[lambda,k, marg(i,7), marg(i,8)]) );
      dzdlambda = ((X.*(lambda*X).^(k-1))/gamma(k)).*exp(-lambda*X) ./ normpdf(zi);
      % Computed dzdk with a forward finite difference scheme
      h = 1/10000;
      cc = ( ferum_cdf(3,X,[lambda,k+h, marg(i,7), marg(i,8)]) - ferum_cdf(3,X,[lambda,k, marg(i,7), marg(i,8)]) ) / h;
      dzdk = cc ./ normpdf(zi);

      dzdmean = dzdk * (2*mean/stdv^2) + dzdlambda * (1/stdv^2) ;
      dzdstdv = dzdk * (-2*(mean^2)/stdv^3) + dzdlambda * (-2*mean/stdv^3) ;

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdlambda;
      J_Z_thetaf.p2    = dzdk;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 4  % Shifted exponential distribution

      mean   = marg(i,2);
      stdv   = marg(i,3);
      lambda = marg(i,5);
      x_zero = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(4,X,[lambda,x_zero, marg(i,7), marg(i,8)]) );
      dzdlambda = (X-x_zero) .* exp(-lambda*(X-x_zero)) ./ normpdf(zi);
      dzdx_zero = -lambda * exp(-lambda*(X-x_zero)) ./ normpdf(zi);

      dzdmean = dzdlambda * 0 + dzdx_zero * 1;
      dzdstdv = dzdlambda * (-1/stdv^2) + dzdx_zero * (-1);

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdlambda;
      J_Z_thetaf.p2    = dzdx_zero;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 5  % Shifted Rayleigh distribution

      mean   = marg(i,2);
      stdv   = marg(i,3);
      a      = marg(i,5);
      x_zero = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(5,X,[a,x_zero, marg(i,7), marg(i,8)]) );
      dzda = -((X-x_zero).^2/a^3).*exp(-0.5*((X-x_zero)/a).^2) ./ normpdf(zi);
      dzdx_zero = -((X-x_zero)/a^2).*exp(-0.5*((X-x_zero)/a).^2) ./ normpdf(zi);

      dzdmean = dzda * 0 + dzdx_zero * 1;
      dzdstdv = dzda * (1/(2-pi/2)^0.5) + dzdx_zero * (-(pi/(4-pi))^0.5);

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzda;
      J_Z_thetaf.p2    = dzdx_zero;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 6  % Uniform distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      a = marg(i,5);
      b = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(6,X,[a,b, marg(i,7), marg(i,8)]) );
      dzda = ((X-b)/(b-a)^2) ./ normpdf(zi);
      dzdb = -((X-a)/(b-a)^2) ./ normpdf(zi);

      dzdmean = dzda * 1 + dzdb * 1;
      dzdstdv = dzda * (-3^0.5) + dzdb * 3^0.5;

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzda;
      J_Z_thetaf.p2    = dzdb;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 7  % Beta distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      q = marg(i,5);
      r = marg(i,6);
      a = marg(i,7);
      b = marg(i,8);

      zi = inv_norm_cdf( ferum_cdf(7,X,[q,r,a,b]) );

      h = 1/10000;

      cca = ( ferum_cdf(7,X,q,r,a+h,b) - ferum_cdf(7,X,[q,r,a,b]) ) / h;
      dzda = cca ./ normpdf(zi);

      ccb = ( ferum_cdf(7,X,q,r,a,b+h) - ferum_cdf(7,X,[q,r,a,b]) ) / h;
      dzdb = ccb ./ normpdf(zi);

      ccq = ( ferum_cdf(7,X,q+h,r,a,b) - ferum_cdf(7,X,[q,r,a,b]) ) / h;
      dzdq = ccq ./ normpdf(zi);

      ccr = ( ferum_cdf(7,X,q,r+h,a,b) - ferum_cdf(7,X,[q,r,a,b]) ) / h;
      dzdr = ccr ./ normpdf(zi);

      dmeandq = (b-a)/(q+r)-q*(b-a)/(q+r)^2;
      dmeandr = -q*(b-a)/(q+r)^2;
      dmeanda = 1-q/(q+r);
      dmeandb = q/(q+r);

      dstdvdq = -(b-a)/(q+r)^2*(q*r/(q+r+1))^(1/2)+1/2*(b-a)/(q+r)/(q*r/(q+r+1))^(1/2)*(r/(q+r+1)-q*r/(q+r+1)^2);
      dstdvdr = -(b-a)/(q+r)^2*(q*r/(q+r+1))^(1/2)+1/2*(b-a)/(q+r)/(q*r/(q+r+1))^(1/2)*(q/(q+r+1)-q*r/(q+r+1)^2);
      dstdvda = -1/(q+r)*(q*r/(q+r+1))^(1/2);
      dstdvdb = 1/(q+r)*(q*r/(q+r+1))^(1/2);

      dzdmean = 0;
      if abs(dmeandq) > eps
         dzdmean = dzdmean + dzdq/dmeandq;
      end
      if abs(dmeandr) > eps
         dzdmean = dzdmean + dzdr/dmeandr;
      end
      if abs(dmeanda) > eps
         dzdmean = dzdmean + dzda/dmeanda;
      end
      if abs(dmeandb) > eps
         dzdmean = dzdmean + dzdb/dmeandb;
      end

      dzdstdv = 0;
      if abs(dstdvdq) > eps
         dzdstdv = dzdstdv + dzdq/dstdvdq;
      end
      if abs(dstdvdr) > eps
         dzdstdv = dzdstdv + dzdr/dstdvdr;
      end
      if abs(dstdvda) > eps
         dzdstdv = dzdstdv + dzda/dstdvda;
      end
      if abs(dstdvdb) > eps
         dzdstdv = dzdstdv + dzdb/dstdvdb;
      end

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdq;
      J_Z_thetaf.p2    = dzdr;
      J_Z_thetaf.p3    = dzda;
      J_Z_thetaf.p4    = dzdb;

   case 8 % Chi-square distribution

      mean   = marg(i,2);
      stdv   = marg(i,3);
      nu     = marg(i,5);
      lambda = 0.5;

      zi = inv_norm_cdf( ferum_cdf(8,X,[nu, marg(i,6), marg(i,7), marg(i,8)]) );
      dzdlambda = 0;
      % Computed dzdnu with a forward finite difference scheme
      h = 1/10000;
      cc = ( ferum_cdf(8,X,[nu+h, marg(i,6), marg(i,7), marg(i,8)]) - ferum_cdf(8,X,[nu, marg(i,6), marg(i,7), marg(i,8)]) ) / h;
      dzdnu = cc ./ normpdf(zi);

      dzdmean = dzdnu * (4*mean/stdv^2);
      dzdstdv = dzdnu * (-4*(mean^2)/stdv^3);

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdnu;
      J_Z_thetaf.p2    = dzdlambda;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 11 % Type I largest value distribution ( same as Gumbel distribution )

      mean = marg(i,2);
      stdv = marg(i,3);
      u_n  = marg(i,5);
      a_n  = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(11,X,[u_n,a_n, marg(i,7), marg(i,8)]) );
      dzdu_n = -a_n * exp(-a_n*(X-u_n)) .* exp(-exp(-a_n*(X-u_n))) ./ normpdf(zi);
      dzda_n = (X-u_n) .* exp(-a_n*(X-u_n)) .* exp(-exp(-a_n*(X-u_n))) ./ normpdf(zi);

      dzdmean = dzdu_n * 1 + dzda_n * 0;
      dzdstdv = dzdu_n * (-0.5772*(6^0.5)/pi) + dzda_n * (-pi/((6^0.5)*stdv^2));

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdu_n;
      J_Z_thetaf.p2    = dzda_n;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 12 % Type I smallest value distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      u_1 = marg(i,5);
      a_1 = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(12,X,[u_1,a_1,marg(i,7), marg(i,8)]) );
      dzdu_1 = -a_1 * exp(a_1*(X-u_1)) .* exp(-exp(a_1*(X-u_1))) ./ normpdf(zi);
      dzda_1 = (X-u_1) .* exp(a_1*(X-u_1)) .* exp(-exp(a_1*(X-u_1))) ./ normpdf(zi);

      dzdmean = dzdu_1 * 1 + dzda_1 * 0;
      dzdstdv = dzdu_1 * 0.5772156649*(6^0.5)/pi + dzda_1 * (-pi/((6^0.5)*stdv^2));

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdu_1;
      J_Z_thetaf.p2    = dzda_1;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 13 % Type II largest value distribution

      mean = marg(i,2);
      stdv = marg(i,3);
      u_n  = marg(i,5);
      k    = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(13,X,[u_n,k, marg(i,7), marg(i,8)]) );
      dzdu_n = (-(u_n./X).^k) * (k/u_n) .* exp(-(u_n./X).^k) ./ normpdf(zi);
      dzdk = (-(u_n./X).^k) .* log(u_n./X) .* exp(-(u_n./X).^k) ./ normpdf(zi);

      J = zeros(2);
      J(1,1) = gamma(1-1/k);
      J(1,2) = (gamma(1-2/k)-(gamma(1-1/k))^2)^0.5;
      h = 1/10000;
      J(2,1) = u_n*( gamma(1-1/(k+h)) - gamma(1-1/k) ) / h;
      J(2,2) = u_n*( (gamma(1-2/(k+h))-(gamma(1-1/(k+h)))^2)^0.5 - (gamma(1-2/k)-(gamma(1-1/k))^2)^0.5 ) / h;
      dzdmeanandstdv = inv(J)*[dzdu_n(:)'; dzdk(:)'];
      dzdmean = reshape(dzdmeanandstdv(1,:),size(X,1),size(X,2));
      dzdstdv = reshape(dzdmeanandstdv(2,:),size(X,1),size(X,2));

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdu_n;
      J_Z_thetaf.p2    = dzdk;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 14 % Type III smallest value distribution

      mean    = marg(i,2);
      stdv    = marg(i,3);
      u_1     = marg(i,5);
      k       = marg(i,6);
      epsilon = marg(i,7);

      zi = inv_norm_cdf( ferum_cdf(14,X,[u_1,k,epsilon, marg(i,8)]) );
      cst1 = ((X - epsilon)/(u_1 - epsilon)).^k ;
      cst2 = exp(-cst1);
      cst3 = (X - epsilon)/(u_1 - epsilon)^2 - 1/(u_1 - epsilon);

      dzdu_1 = -cst1 .* cst2 * (k/(u_1-epsilon)) ./ normpdf(zi);
      dzdk = cst1 .* cst2 .* log((X - epsilon)/(u_1 - epsilon)) ./ normpdf(zi);
      dzdepsilon = cst1 .* cst2 * k .* cst3 * (u_1 - epsilon) ./ ( (X-epsilon) .* normpdf(zi) );

      dmeandu_1 = gamma(1+1/k);
      dmeandk = -(u_1-epsilon)*Psi(1+1/k)*gamma(1+1/k)/k^2;
      dmeandepsilon = 1-gamma(1+1/k);

      dstdvdu_1 = (gamma(1+2/k)-gamma(1+1/k)^2)^(1/2);
      dstdvdk = 1/2*(u_1-epsilon)/(gamma(1+2/k)-gamma(1+1/k)^2)^(1/2)*(-2*Psi(1+2/k)*gamma(1+2/k)/k^2+2*gamma(1+1/k)^2*Psi(1+1/k)/k^2);
      dstdvdepsilon = -(gamma(1+2/k)-gamma(1+1/k)^2)^(1/2);

      dzdmean = 0;
      if abs(dmeandu_1) > eps
         dzdmean = dzdmean + dzdu_1/dmeandu_1;
      end
      if abs(dmeandk) > eps
         dzdmean = dzdmean + dzdk/dmeandk;
      end
      if abs(dmeandepsilon) > eps
         dzdmean = dzdmean + dzdepsilon/dmeandepsilon;
      end

      dzdstdv = 0;
      if abs(dstdvdu_1) > eps
         dzdstdv = dzdstdv + dzdu/dstdvdu_1;
      end
      if abs(dstdvdk) > eps
         dzdstdv = dzdstdv + dzdk/dstdvdk;
      end
      if abs(dstdvdepsilon) > eps
         dzdstdv = dzdstdv + dzdepsilon/dstdvdepsilon;
      end

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdu_1;
      J_Z_thetaf.p2    = dzdk;
      J_Z_thetaf.p3    = dzdepsilon;
      J_Z_thetaf.p4    = 0;

   case 15 % Gumbel distribution ( same as type I largest value distribution )

      mean = marg(i,2);
      stdv = marg(i,3);
      u_n  = marg(i,5);
      a_n  = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(15,X,[u_n,a_n, marg(i,7), marg(i,8)]) );
      dzdu_n = -a_n * exp(-a_n*(X-u_n)) .* exp(-exp(-a_n*(X-u_n))) ./ normpdf(zi);
      dzda_n = (X-u_n) .* exp(-a_n*(X-u_n)) .* exp(-exp(-a_n*(X-u_n))) ./ normpdf(zi);

      dzdmean = dzdu_n * 1 + dzda_n * 0;
      dzdstdv = dzdu_n * (-0.5772*(6^0.5)/pi) + dzda_n * (-pi/((6^0.5)*stdv^2));

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdu_n;
      J_Z_thetaf.p2    = dzda_n;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )

      mean = marg(i,2);
      stdv = marg(i,3);
      u_1  = marg(i,5);
      k    = marg(i,6);

      zi = inv_norm_cdf( ferum_cdf(16,X,[u_1,k, marg(i,7), marg(i,8)]) );
      dzdu_1 = -((X/u_1).^k) * (k/u_1) .* exp(-(X/u_1).^k) ./ normpdf(zi);
      dzdk = ((X/u_1).^k) .* log(X/u_1) .* exp(-(X/u_1).^k) ./ normpdf(zi);

      J = zeros(2);
      J(1,1) = gamma(1+1/k);
      J(1,2) = (gamma(1+2/k)-(gamma(1+1/k))^2)^0.5;
      h = 1/10000;
      J(2,1) = u_1*( gamma(1+1/(k+h)) - gamma(1+1/k) ) / h;
      J(2,2) = u_1*( (gamma(1+2/(k+h))-(gamma(1+1/(k+h)))^2)^0.5 - (gamma(1+2/k)-(gamma(1+1/k))^2)^0.5 ) / h;
      dzdmeanandstdv = inv(J)*[dzdu_1(:)'; dzdk(:)'];
      dzdmean = reshape(dzdmeanandstdv(1,:),size(X,1),size(X,2));
      dzdstdv = reshape(dzdmeanandstdv(2,:),size(X,1),size(X,2));

      J_Z_thetaf.mu    = dzdmean;
      J_Z_thetaf.sigma = dzdstdv;
      J_Z_thetaf.p1    = dzdu_1;
      J_Z_thetaf.p2    = dzdk;
      J_Z_thetaf.p3    = 0;
      J_Z_thetaf.p4    = 0;

   case 18 % Laplace distribution

   case 19 % Pareto distribution

   case 51 % Truncated normal marginal distribution

      mean = marg(i,5);
      stdv = marg(i,6);
      xmin = marg(i,7);
      xmax = marg(i,8);

      zi = inv_norm_cdf( ferum_cdf(51,X,[mean,stdv,xmin,xmax]) );

      h = 1/10000;

      ccmean = ( ferum_cdf(51,X,[mean+h,stdv,xmin,xmax]) - ferum_cdf(51,X,[mean,stdv,xmin,xmax]) ) / h;
      dzdmean = ccmean ./ normpdf(zi);
      ccstdv = ( ferum_cdf(51,X,[mean,stdv+h,xmin,xmax]) - ferum_cdf(51,X,[mean,stdv,xmin,xmax]) ) / h;
      dzdstdv = ccstdv ./ normpdf(zi);
      ccxmin = ( ferum_cdf(51,X,[mean,stdv,xmin+h,xmax]) - ferum_cdf(51,X,[mean,stdv,xmin,xmax]) ) / h;
      dzdxmin = ccxmin ./ normpdf(zi);
      ccxmax = ( ferum_cdf(51,X,[mean,stdv,xmin,xmax]) - ferum_cdf(51,X,[mean,stdv,xmin,xmax-h]) ) / h;
      dzdxmax = ccxmax ./ normpdf(zi);

      h = 1/10000;

      dmean_tdmean = ( mean_norm_truncated([mean+h,stdv,xmin,xmax]) - mean_norm_truncated([mean,stdv,xmin,xmax]) ) / h;
      dmean_tdstdv = ( mean_norm_truncated([mean,stdv+h,xmin,xmax]) - mean_norm_truncated([mean,stdv,xmin,xmax]) ) / h;
      dmean_tdxmin = ( mean_norm_truncated([mean,stdv,xmin+h,xmax]) - mean_norm_truncated([mean,stdv,xmin,xmax]) ) / h;
      dmean_tdxmax = ( mean_norm_truncated([mean,stdv,xmin,xmax+h]) - mean_norm_truncated([mean,stdv,xmin,xmax]) ) / h;

      dstdv_tdmean = ( stdv_norm_truncated([mean+h,stdv,xmin,xmax]) - stdv_norm_truncated([mean,stdv,xmin,xmax]) ) / h;
      dstdv_tdstdv = ( stdv_norm_truncated([mean,stdv+h,xmin,xmax]) - stdv_norm_truncated([mean,stdv,xmin,xmax]) ) / h;
      dstdv_tdxmin = ( stdv_norm_truncated([mean,stdv,xmin+h,xmax]) - stdv_norm_truncated([mean,stdv,xmin,xmax]) ) / h;
      dstdv_tdxmax = ( stdv_norm_truncated([mean,stdv,xmin,xmax+h]) - stdv_norm_truncated([mean,stdv,xmin,xmax]) ) / h;

      J_Z_thetaf.mu = 0;
      if abs(dmean_tdmean) > eps
         J_Z_thetaf.mu = J_Z_thetaf.mu + dzdmean/dmean_tdmean;
      end
      if abs(dmean_tdstdv) > eps
         J_Z_thetaf.mu = J_Z_thetaf.mu + dzdstdv/dmean_tdstdv;
      end
      if abs(dmean_tdxmin) > eps
         J_Z_thetaf.mu = J_Z_thetaf.mu + dzdxmin/dmean_tdxmin;
      end
      if abs(dmean_tdxmax) > eps
         J_Z_thetaf.mu = J_Z_thetaf.mu + dzdxmax/dmean_tdxmax;
      end

      J_Z_thetaf.sigma = 0;
      if abs(dstdv_tdmean) > eps
         J_Z_thetaf.sigma = J_Z_thetaf.sigma + dzdmean/dstdv_tdmean;
      end
      if abs(dstdv_tdstdv) > eps
         J_Z_thetaf.sigma = J_Z_thetaf.sigma + dzdstdv/dstdv_tdstdv;
      end
      if abs(dstdv_tdxmin) > eps
         J_Z_thetaf.sigma = J_Z_thetaf.sigma + dzdxmin/dstdv_tdxmin;
      end
      if abs(dstdv_tdxmax) > eps
         J_Z_thetaf.sigma = J_Z_thetaf.sigma + dzdxmax/dstdv_tdxmax;
      end

      J_Z_thetaf.p1 = dzdmean;
      J_Z_thetaf.p2 = dzdstdv;
      J_Z_thetaf.p3 = dzdxmin;
      J_Z_thetaf.p4 = dzdxmax;

   otherwise

end