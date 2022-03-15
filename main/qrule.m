function [bp,wf]=qrule(n,wfun,alpha,beta)
%QRULE compute abscissas and weight factors for Gaussian quadratures
%
%CALL:  [bp,wf]=qrule(n,wfun,alpha,beta)
%  
%  bp = base points (abscissas)
%  wf = weight factors
%  n  = number of base points (abscissas) (integrates a (2n-1)th order
%       polynomial exactly)
%wfun = weight function%     
%     1  p(x)=1                       a =-1,   b = 1 Legendre (default)
%     2  p(x)=1/sqrt((x-a)*(b-x)),    a =-1,   b = 1 Chebyshev of the
%                                                             first kind
%     3  p(x)=sqrt((x-a)*(b-x)),      a =-1,   b = 1 Chebyshev of the 
%                                                            second kind
%     4  p(x)=sqrt((x-a)/(b-x)),      a = 0,   b = 1
%     5  p(x)=1/sqrt(b-x),            a = 0,   b = 1
%     6  p(x)=sqrt(b-x),              a = 0,   b = 1
%     7  p(x)=(x-a)^alpha*(b-x)^beta  a =-1,   b = 1 Jacobi 
%                                     alpha, beta >-1 (default alpha=beta=0)
%     8  p(x)=x^alpha*exp(-x)         a = 0,   b = inf generalized Laguerre
%     9  p(x)=exp(-x^2)               a =-inf, b = inf Hermite
%    10  p(x)=1                       a =-1,   b = 1 Legendre (slower than 1)
%
%  The Gaussian Quadrature integrates a (2n-1)th order
%  polynomial exactly and the integral is of the form
%           b                         n
%          Int ( p(x)* F(x) ) dx  =  Sum ( wf_j* F( bp_j ) )
%           a                        j=1		          
%  See also: gaussq

% Reference 
%   wfun 1: copied from grule.m in NIT toolbox, see ref [2] 
%   wfun 2-6: see ref [4]
%   wfun 7-10:  Adapted from Netlib routine gaussq.f see ref [1,3]
%
% [1]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% [4] Abromowitz and Stegun (1954) ''

%  By Bryce Gardner, Purdue University, Spring 1993.
% Modified by Per A. Brodtkorb 19.02.99 pab@marin.ntnu.no
% to compute other quadratures  than the default
if nargin<4|isempty(beta),
 beta=0; 
end

if nargin<3|isempty(alpha),
  alpha=0; 
end
if alpha<=-1 | beta <=-1,
  error('alpha and beta must be greater than -1')
end

if nargin<2|isempty(wfun),
  wfun=1; 
end	


switch wfun, %
  case 1,
    %  This routine computes Gauss base points and weight factors
    %  using the algorithm given by Davis and Rabinowitz in 'Methods
    %  of Numerical Integration', page 365, Academic Press, 1975.
    bp=zeros(n,1); wf=bp; iter=2; m=fix((n+1)/2); e1=n*(n+1);
    mm=4*m-1; t=(pi/(4*n+2))*(3:4:mm); nn=(1-(1-1/n)/(8*n*n));
    xo=nn*cos(t);
    for j=1:iter
      pkm1=1; pk=xo;
      for k=2:n
	t1=xo.*pk; pkp1=t1-pkm1-(t1-pkm1)/k+t1;
	pkm1=pk; pk=pkp1;
      end
      den=1.-xo.*xo; d1=n*(pkm1-xo.*pk); dpn=d1./den;
      d2pn=(2.*xo.*dpn-e1.*pk)./den;
      d3pn=(4*xo.*d2pn+(2-e1).*dpn)./den;
      d4pn=(6*xo.*d3pn+(6-e1).*d2pn)./den;
      u=pk./dpn; v=d2pn./dpn;
      h=-u.*(1+(.5*u).*(v+u.*(v.*v-u.*d3pn./(3*dpn))));
      p=pk+h.*(dpn+(.5*h).*(d2pn+(h/3).*(d3pn+.25*h.*d4pn)));
      dp=dpn+h.*(d2pn+(.5*h).*(d3pn+h.*d4pn/3));
      h=h-p./dp; xo=xo+h;
    end
    bp=-xo-h;
    fx=d1-h.*e1.*(pk+(h/2).*(dpn+(h/3).*(...
	d2pn+(h/4).*(d3pn+(.2*h).*d4pn))));
    wf=2*(1-bp.^2)./(fx.*fx);
    if (m+m) > n, bp(m)=0; end
    if ~((m+m) == n), m=m-1; end
    jj=1:m; n1j=(n+1-jj); bp(n1j)=-bp(jj); wf(n1j)=wf(jj);
    % end
    
 case 2, % p(x)=1/sqrt((x-a)*(b-x)), a=-1 and b=1 (default) 
  j=[1:n];
  wf = ones(1,n) * pi / n;
  bp=cos( (2*j-1)*pi / (2*n) );

 case 3, %p(x)=sqrt((x-a)*(b-x)),   a=-1   and b=1
  j=[1:n];
  wf = pi/ (n+1) *sin( j*pi / (n+1) ).^2;
  bp=cos( j*pi / (n+1) );

 case 4, %p(x)=sqrt((x-a)/(b-x)),   a=0   and b=1
    j=[1:n];
    bp=cos( (2*j-1)*pi /2/ (2*n+1) ).^2;
    wf=2*pi.*bp/(2*n+1) ;

 case 5, % %p(x)=1/sqrt(b-x),   a=0   and b=1
   [bp wf]=grule(2*n);
  wf(bp<0)=[];
  wf=wf*2;
   bp(bp<0)=[];
  bp=1-bp.^2;

 case 6, % %p(x)=sqrt(b-x),   a=0   and b=1
   [bp wf]=grule(2*n+1);
  wf(bp<=0)=[];
   bp(bp<=0)=[];
  wf=2*bp.^2.*wf;
  bp=1-bp.^2;
  
 case {7,8,9,10} ,%
  %7 p(x)=(x-a)^alpha*(b-x)^beta a=-1 b=1 Jacobi
  %8 p(x)=x^alpha*exp(-x) a=0,   b=inf generalized Laguerre
  %9 p(x)=exp(-x^2)       a=-inf, b=inf Hermite 
  %10 p(x)=1               a=-1 b=1        Legendre slower than 1
  % this procedure uses the coefficients a(j), b(j) of the
  %      recurrence relation
  %
  %           b p (x) = (x - a ) p   (x) - b   p   (x)
  %            j j            j   j-1       j-1 j-2
  %
  %      for the various classical (normalized) orthogonal polynomials,
  %      and the zero-th moment
  %
  %           muzero = integral w(x) dx
  %
  %      of the given polynomial's weight function w(x).  since the
  %      polynomials are orthonormalized, the tridiagonal matrix is
  %      guaranteed to be symmetric.
  %
  % 
  %         the input parameter alpha is used only for laguerre and
  %      jacobi polynomials, and the parameter beta is used only for
  %      jacobi polynomials.  the laguerre and jacobi polynomials
  %      require the gamma function.

  a=zeros(n,1);
  b=zeros(n-1,1);
  switch wfun
    case 7,  %jacobi
      ab = alpha + beta;
      abi = 2 + ab;
      muzero = 2^(ab + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(abi);
      a(1) = (beta - alpha)/abi;
      b(1) = sqrt(4*(1 + alpha)*(1 + beta)/((abi + 1)*abi^2));
      a2b2 = beta^2 - alpha^2;
      
      i = (2:n-1)';
      abi = 2*i + ab;
      a(i) = a2b2./((abi - 2).*abi);
      a(n) =a2b2./((2*n - 2+ab).*(2*n+ab));
      b(i) = sqrt (4*i.*(i + alpha).*(i + beta)*(i + ab)./((abi.^2 - 1).*abi.^2));
   
    case 8, % Laguerre
      muzero=gamma(alpha+1);
      i = (1:n-1)';
      a(i) = 2 .* i - 1 + alpha;
      a(n)=2*n-1+alpha;
      b = sqrt( i .* (i + alpha) );
    case 9, %Hermite 
      i = (1:(n-1))';
      muzero = sqrt(pi);
      %a=zeros(m,1);
      b=sqrt(i/2);    
    case 10,  % legendre NB! much slower than wfun=1        
      muzero = 2;
      i = (1:n-1)';
      abi = i;
      b(i) = abi./sqrt(4*abi.^2 - 1);
      
  end
   
  %[v d] = eig( full(spdiags([b a b],-1:1,n,n )));
  [v d ] = eig( diag(a) + diag(b,1) + diag(b,-1) );
  wf = v(1,:);
  if 1,
    [bp i] = sort( diag(d) );
    wf = wf(i);
  else % save some valuable time by not sorting
    bp = diag(d) ;
  end
  bp=bp';
  
  wf = muzero.* wf.^2;

otherwise, error('unknown weight function')
end

% end
