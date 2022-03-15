function  J_u_x = jacobian(x,u,probdata)

transf_type = probdata.transf_type;
marg = probdata.marg;
nrv = size(marg,1);

%! NOTE!
%! the same procedure is done in every cases, this whole case stuff could be replaced
% exp(log(pdf1) - log(pdf2)) might be numerically better, log pdf could be
% an option in the functions

switch transf_type
    
    case 1
        
    case 2
        
        %J_u_x = zeros(nrv);
        %
        %for i = 1 : nrv
        %   switch marg(i,1)
        %      case 1   % Normal distribution
        %         J_u_x(i,i) = 1/marg(i,3);
        %      case 2   % Lognormal distribution
        %         ksi = sqrt( log( 1 + ( marg(i,3) / marg(i,2) )^2 ) );
        %         J_u_x(i,i) = 1 / ( ksi * x(i) );
        %      case 4  % Shifted exponential marginal distribution
        %         pdf1 = ferum_pdf(4,x(i),marg(i,2),marg(i,3));
        %         pdf2 = ferum_pdf(1,u(i),0,1);
        %         J_u_x(i,i) = pdf1/pdf2;
        %      case 5  % Shifted Rayleigh marginal distribution
        %         pdf1 = ferum_pdf(5,x(i),marg(i,2),marg(i,3));
        %         pdf2 = ferum_pdf(1,u(i),0,1);
        %         J_u_x(i,i) = pdf1/pdf2;
        %      case 6  % Uniform marginal distribution
        %         pdf1 = ferum_pdf(6,x(i),marg(i,2),marg(i,3));
        %         pdf2 = ferum_pdf(1,u(i),0,1);
        %         J_u_x(i,i) = pdf1/pdf2;
        %      case 11  % Type I Largest Value or Gumbel marginal distribution
        %         pdf1 = ferum_pdf(11,x(i),marg(i,2),marg(i,3));
        %         pdf2 = ferum_pdf(1,u(i),0,1);
        %         J_u_x(i,i) = pdf1/pdf2;
        %      otherwise
        %   end
        %end
        
    case 3
        
        Lo = probdata.Lo;
        iLo = probdata.iLo;
        z = Lo * u;
        J_z_x = zeros(nrv);
        
        for i = 1 : nrv
            switch marg(i,1)
                case 1   % Normal distribution
                    pdf1 = ferum_pdf(1,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
%                     pdf2 = normpdf_mp(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                    %J_z_x(i,i) = 1/marg(i,3);
                case 2   % Lognormal distribution
                    pdf1 = ferum_pdf(2,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                    %                ksi = sqrt( log( 1 + ( marg(i,3) / marg(i,2) )^2 ) );
                    %                J_z_x(i,i) = 1 / ( ksi * x(i) );
                case 3  % Gamma distribution
                    %lambda = marg(i,5);
                    %k = marg(i,6);
                    %pdf1 = lambda * (lambda*x(i))^(k-1) / gamma(k) * exp(-lambda*x(i));
                    pdf1 = ferum_pdf(3,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 4  % Shifted exponential distribution
                    %lambda = marg(i,5);
                    %x_zero = marg(i,6);
                    %pdf1 = lambda * exp(-lambda*(x(i)-x_zero));
                    pdf1 = ferum_pdf(4,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 5  % Shifted Rayleigh distribution
                    %a = marg(i,5);
                    %x_zero = marg(i,6);
                    %pdf1 = (x(i)-x_zero)/a^2 * exp(-0.5*((x(i)-x_zero)/a)^2);
                    pdf1 = ferum_pdf(5,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 6  % Uniform distribution
                    %a = marg(i,5);
                    %b = marg(i,6);
                    %pdf1 = 1 / (b-a);
                    pdf1 = ferum_pdf(6,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 7  % Beta distribution
                    %q = marg(i,5);
                    %r = marg(i,6);
                    %a = marg(i,7);
                    %b = marg(i,8);
                    %pdf1 = (x(i)-a)^(q-1) * (b-x(i))^(r-1) / ( (gamma(q)*gamma(r)/gamma(q+r)) * (b-a)^(q+r-1) );
                    pdf1 = ferum_pdf(7,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 8 % Chi-square distribution
                    %lambda = 0.5;
                    %nu = marg(i,5);
                    %k = nu/2;
                    %pdf1 = lambda * (lambda*x(i))^(k-1) * exp(-lambda*x(i)) / gamma(k) ;
                    pdf1 = ferum_pdf(8,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i)= pdf1/pdf2;
                case 11 % Type I largest value distribution ( same as Gumbel distribution )
                    %u_n = marg(i,5);
                    %a_n = marg(i,6);
                    %pdf1 = a_n * exp( -a_n*(x(i)-u_n) - exp(-a_n*(x(i)-u_n)) );
                    pdf1 = ferum_pdf(11,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
%                     pdf2 = normpdf_mp(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 12 % % Type I smallest value distribution
                    %u_1 = marg(i,5);
                    %a_1 = marg(i,6);
                    %pdf1 = a_1 * exp( a_1*(x(i)-u_1) - exp(a_1*(x(i)-u_1)) );
                    pdf1 = ferum_pdf(12,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 13 % Type II largest value distribution
                    %u_n = marg(i,5);
                    %k = marg(i,6);
                    %pdf1 = k/u_n * (u_n/x(i))^(k+1) * exp(-(u_n/x(i))^k);
                    pdf1 = ferum_pdf(13,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 14 % Type III smallest value distribution
                    %u_1 = marg(i,5);
                    %k = marg(i,6);
                    %epsilon = marg(i,7);
                    %pdf1 = k/(u_1-epsilon) * ((x(i)-epsilon)/(u_1-epsilon))^(k-1) ...
                    %       * exp(-((x(i)-epsilon)/(u_1-epsilon))^k);
                    pdf1 = ferum_pdf(14,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 15 % Gumbel distribution ( same as type I largest value distribution )
                    %u_n = marg(i,5);
                    %a_n = marg(i,6);
                    %pdf1 = a_n * exp( -a_n*(x(i)-u_n) - exp(-a_n*(x(i)-u_n)) );
                    pdf1 = ferum_pdf(15,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 16 % Weibull distribution ( same as Type III smallest value distribution with epsilon = 0 )
                    %u_1 = marg(i,5);
                    %k = marg(i,6);
                    %pdf1 = k/u_1 * (x(i)/u_1)^(k-1) * exp(-(x(i)/u_1)^k);
                    pdf1 = ferum_pdf(16,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 20
                    pdf1 = ferum_pdf(20, x(i), marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 25
                    pdf1 = ferum_pdf(25, x(i), marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 30
                    pdf1 = ferum_pdf(30, x(i), marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                case 31
                    pdf1 = ferum_pdf(31, x(i), marg(i,5:8));
                    pdf2 = normpdf(z(i));
%                     if isnan(pdf1/pdf2)
%                         keyboard
%                     end
                    J_z_x(i,i) = pdf1/pdf2;
                case 51  % Truncated normal marginal distribution
                    pdf1 = ferum_pdf(51,x(i),marg(i,5:8));
                    pdf2 = normpdf(z(i));
                    J_z_x(i,i) = pdf1/pdf2;
                otherwise
            end
        end
        
        J_u_x = iLo * J_z_x;
        
    case 4
        
        %alpha12 = 1;
        %
        %J_u_x = zeros(2);
        %
        %p1 = ferum_pdf(marg(1,1),x(1),marg(1,2),marg(1,3));
        %p2 = ferum_pdf(marg(2,1),x(2),marg(2,2),marg(2,3));
        %
        %P1 = ferum_cdf(marg(1,1),x(1),marg(1,2),marg(1,3));
        %P2 = ferum_cdf(marg(2,1),x(2),marg(2,2),marg(2,3));
        %
        %J_u_x(1,1) = pdf_morgenstern(x,marg,1) / ferum_pdf(1,u(1),0,1);
        %J_u_x(2,2) = pdf_morgenstern(x,marg,2) / ferum_pdf(1,u(2),0,1);
        %
        %J_u_x(2,1) = -2 * alpha12 * P2 * (1-P2) * p1 / ferum_pdf(1,u(2),0,1);
        
    otherwise
        
end