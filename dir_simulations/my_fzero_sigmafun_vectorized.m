function [b,fval,exitflag,output] = my_fzero_sigmafun_vectorized(x,tolx,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield)

% fzero search around X = <x1 ... xn>
% allr = my_fzero_sigmafun_vectorized(X,tolx,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
%
% fzero search in [a b] interval
% allr = my_fzero_sigmafun_vectorized([a b]'*ones(1,n),tolx,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);



%function [b,fval,exitflag,output] = my_fzero(FunFcnIn,x,options,varargin)
%
%FZERO  Scalar nonlinear zero finding. 
%   X = FZERO(FUN,X0) tries to find a zero of the function FUN near X0, 
%   if X0 is a scalar.  It first finds an interval containing X0 where the 
%   function values of the interval endpoints differ in sign, then searches 
%   that interval for a zero. FUN accepts real scalar input X and returns 
%   a real scalar function value F, evaluated at X.  The value X returned 
%   by FZERO is near a point where FUN changes sign (if FUN is continuous), 
%   or NaN if the search fails.  
%
%   X = FZERO(FUN,X0), where X0 is a vector of length 2, assumes X0 is an 
%   interval where the sign of FUN(X0(1)) differs from the sign of FUN(X0(2)).
%   An error occurs if this is not true.  Calling FZERO with an interval  
%   guarantees FZERO will return a value near a point where FUN changes 
%   sign.
%
%   X = FZERO(FUN,X0), where X0 is a scalar value, uses X0 as a starting 
%   guess. FZERO looks for an interval containing a sign change for FUN and 
%   containing X0.  If no such interval is found, NaN is returned.  
%   In this case, the search terminates when the search interval 
%   is expanded until an Inf, NaN, or complex value is found. Note: if
%   the option FunValCheck is 'on', then an error will occur if an NaN or 
%   complex value is found.
%
%   X = FZERO(FUN,X0,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details.  Used
%   options are Display, TolX, FunValCheck, and OutputFcn. Use OPTIONS = [] 
%   as a place holder if no options are set.
%
%   [X,FVAL]= FZERO(FUN,...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FZERO(...) returns an EXITFLAG that describes the 
%   exit condition of FZERO. Possible values of EXITFLAG and the corresponding 
%   exit conditions are
%
%     1  FZERO found a zero X.
%    -1  Algorithm terminated by output function.
%    -3  NaN or Inf function value encountered during search for an interval
%         containing a sign change.
%    -4  Complex function value encountered during search for an interval 
%         containing a sign change.
%    -5  FZERO may have converged to a singular point.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FZERO(...) returns a structure OUTPUT
%   with the number of function evaluations in OUTPUT.funcCount, the
%   algorithm name in OUTPUT.algorithm, the number of iterations to
%   find an interval (if needed) in OUTPUT.intervaliterations, the
%   number of zero-finding iterations in OUTPUT.iterations, and the
%   exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fzero(@sin,3)
%     returns pi.
%        X = fzero(@sin,3,optimset('disp','iter')) 
%     returns pi, uses the default tolerance and displays iteration information.
%
%     FUN can also be an anonymous function:
%        X = fzero(@(x) sin(3*x),2)
%
%   If FUN is parameterized, you can use anonymous functions to capture the 
%   problem-dependent parameters. Suppose you want to solve the equation given
%   in the function MYFUN, which is parameterized by its second argument A. Here
%   MYFUN is an M-file function such as
%
%       function f = myfun(x,a)
%       f = cos(a*x);
%
%   To solve the equation for a specific value of A, first assign the value to A.
%   Then create a one-argument anonymous function that captures that value of A 
%   and calls MYFUN with two arguments. Finally, pass this anonymous function to 
%   FZERO:
%
%       a = 2; % define parameter first
%       x = fzero(@(x) myfun(x,a),0.1)
%   
%   Limitations
%        X = fzero(@(x) abs(x)+1, 1) 
%     returns NaN since this function does not change sign anywhere on the 
%     real axis (and does not have a zero as well).
%        X = fzero(@tan,2)
%     returns X near 1.5708 because the discontinuity of this function near the 
%     point X gives the appearance (numerically) that the function changes sign at X.
%
%   See also ROOTS, FMINBND, FUNCTION_HANDLE.

%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.33.4.7 $  $Date: 2004/04/16 22:05:26 $

%  This algorithm was originated by T. Dekker.  An Algol 60 version,
%  with some improvements, is given by Richard Brent in "Algorithms for
%  Minimization Without Derivatives", Prentice-Hall, 1973.  A Fortran
%  version is in Forsythe, Malcolm and Moler, "Computer Methods
%  for Mathematical Computations", Prentice-Hall, 1976.

global nfun


% Initialization
fcount = 0;
iter = 0;
intervaliter = 0;
exitflag = 1;
procedure = ' ';

nx = size(alla,2);

tol = tolx;


if isfield(analysisopt,'sigmafun_write2file')
   if analysisopt.sigmafun_write2file
      fid = fopen('my_fzero_sigmafun.txt','wt');
      printtype = 'iter';
   else
      fid = 1;
      printtype = 'off';
   end
else
   fid = 1;
   printtype = 'off';
end
switch printtype
    case 'notify'
        trace = 1;
    case {'none', 'off'}
        trace = 0;
    case 'iter'
        trace = 3;
    case 'debug'
        trace = 4;
    case 'final'
        trace = 2;
    otherwise
        trace = 1;
end

% Interval input
if (size(x,1) == 2)

    if trace > 2
        fprintf(fid,'\n'); %Initial blank line
    end

    a = x(1,:); savea=a;
    b = x(2,:); saveb=b;

    fa = sigmafun(a,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
    fb = sigmafun(b,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
    fcount = fcount + 2*nx;

    if any(~isfinite([fa fb])) || any(~isreal([fa fb]))
        error('MATLAB:fzero:ValuesAtEndPtsComplexOrNotFinite',...
            'Function values at interval endpoints must be finite and real.')
    end

    savefa = fa; savefb = fb;
    
    if ( fa == 0 )
        b = a;
        msg = sprintf('Zero find terminated.\n');
        if trace > 1
            fprintf(fid,msg);
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fa;
        return
    elseif ( fb == 0)
        % b = b;
        msg = sprintf('Zero find terminated.\n');
        if trace > 1
            fprintf(fid,msg);
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fb;
        return
    elseif (fa > 0) == (fb > 0)
        error('MATLAB:fzero:ValuesAtEndPtsSameSign',...
            'The function values at the interval endpoints must differ in sign.')
    end
    
    % Starting guess scalar input

elseif (size(x,1) == 1)

    if trace > 2 
        fprintf(fid,'\n');
        fprintf(fid,'Search for an interval around %.4e containing a sign change:\n', x);
        header = '         left       f(left)         right      f(right)\n';
    end

    fx = sigmafun(x,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
    fcount = fcount + 1*nx;
    
    if fx == 0
        b = x;
        msg = sprintf('Zero find terminated.\n');
        if trace > 1
            fprintf(fid,msg);
        end
        output.intervaliterations = intervaliter;
        output.iterations = iter;
        output.funcCount = fcount;
        output.algorithm = 'bisection, interpolation';
        output.message = msg;
        fval = fx;
        return
    elseif any(~isfinite(fx)) | any(~isreal(fx))
        error('MATLAB:fzero:ValueAtInitGuessComplexOrNotFinite',...
            'Function value at starting guess must be finite and real.');
    end

    dx = x/50;
    i_zeros = find(x==0);
    dx(i_zeros) = 1/50;
    
    % Find change of sign.
    twosqrt = sqrt(2)*ones(1,nx); 
    a = x; fa = fx; b = x; fb = fx;
    
    if trace > 2
        procedure='initial interval';
        fprintf(fid,'Func-count = %5.0f   -   Procedure : %s\n', fcount, procedure);
        fprintf(fid,header);
        fprintf(fid,'%13.6g %13.6g %13.6g %13.6g\n',[a;fa;b;fb]);
    end

    while any ( (fa > 0) == (fb > 0) )
        
        i_diff_sign = find( ~( (fa > 0) == (fb > 0) ) );
        if ~isempty(i_diff_sign)
           twosqrt(i_diff_sign) = 1;
        end
       
        intervaliter = intervaliter + 1;

        dx = twosqrt.*dx;
        
        a = x - dx;  fa = sigmafun(a,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
        fcount = fcount + 1*nx;
        if any(~isfinite(fa)) | any(~isreal(fa))
            [exitflag,msg] = disperr(a,fa,trace);
            b = NaN; fval = NaN;
            output.intervaliter = intervaliter;
            output.iterations = iter;
            output.funcCount = fcount;
            output.algorithm = 'bisection, interpolation';
            output.message = msg;
            return
        end

        if all ( (fa > 0) ~= (fb > 0) ) % check for different sign
            % Before we exit the while loop, print out the latest interval
            if trace > 2
                procedure='search';
                fprintf(fid,'Func-count = %5.0f   -   Procedure : %s\n', fcount, procedure);
                fprintf(fid,header);
                fprintf(fid,'%13.6g %13.6g %13.6g %13.6g\n',[a;fa;b;fb]);
            end
            break
        end
        
        b = x + dx;  fb = sigmafun(b,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
        fcount = fcount + 1*nx;        
        if any(~isfinite(fb)) | any(~isreal(fb))
            [exitflag,msg] = disperr(b,fb,trace);
            b = NaN; fval = NaN;
            output.intervaliter = intervaliter;
            output.iterations = iter;
            output.funcCount = fcount;
            output.algorithm = 'bisection, interpolation';
            output.message = msg;
            return
        end

        if trace > 2
            procedure='search';
            fprintf(fid,'Func-count = %5.0f   -   Procedure : %s\n', fcount, procedure);
            fprintf(fid,header);
            fprintf(fid,'%13.6g %13.6g %13.6g %13.6g\n',[a;fa;b;fb]);
        end

    end % while
    
    if trace > 2
        fprintf(fid,'\n');
        fprintf(fid,'Search for a zero in the interval [%13.6g, %13.6g]\n',[a;b]);
    end
    
    savea = a; savefa = fa; saveb = b; savefb = fb;

else

    error('MATLAB:fzero:LengthArg2', 'Second argument must be of length 1 or 2.');

end % if (length(x) == 2



c = b; % Added by JMB
fc = fb;

d = b - a; % Added by JMB
e = d; % Added by JMB

%procedure = 'initial';
procedure = zeros(1,nx);
header2 = '            x          f(x)     Procedure (0: init, 1:bisect, 2:interp)\n';

% Main loop, exit from middle of the loop
while any(fb ~= 0)

    % Insure that b is the best result so far, a is the previous
    % value of b, and c is on the opposite of the zero from b.

    i_c = find( (fb > 0) == (fc > 0) );
    if ~isempty(i_c)
       c(i_c) = a(i_c); fc(i_c) = fa(i_c);
       d(i_c) = b(i_c) - a(i_c);  e(i_c) = d(i_c);
       if trace > 3
           fprintf(fid,'#01# c = %13.6g  -  fc = %13.6g  -  d = %13.6g  -  e = %13.6g\n', [c; fc; d; e]);
       end
    end

    i_abc = find( abs(fc) < abs(fb) );
    if ~isempty(i_abc)
       a(i_abc) = b(i_abc);    b(i_abc) = c(i_abc);    c(i_abc) = a(i_abc);
       fa(i_abc) = fb(i_abc);  fb(i_abc) = fc(i_abc);  fc(i_abc) = fa(i_abc);
       if trace > 3
          fprintf(fid,'#02# a = %13.6g  -  fa = %13.6g  -  b = %13.6g  -  fb = %13.6g  -  c = %13.6g  -  fc = %13.6g\n', [a; fa; b; fb; c; fc]);
      end
    end
    
    % Convergence test and possible exit
    m = 0.5*(c - b);
    toler = 2.0*tol*max(abs(b),1.0);
    if trace > 3
        fprintf(fid,'#03# m = %13.6g  -  toler = %13.6g\n', [m; toler]);
    end
    if all( (abs(m) <= toler) ) | all( fb == 0.0 )
        if trace > 3
            fprintf(fid,'#04# all( (abs(m) <= toler) ) -> %d  -  all( fb == 0.0 ) -> %d\n', all( (abs(m) <= toler) ), all( fb == 0.0 ));
        end
        break
    end

    if trace > 2
        fprintf(fid,'Func-count = %5.0f\n', fcount);
        fprintf(fid,header2);
        fprintf(fid,'%13.6g %13.6g %13.6g\n',[b;fb;procedure]);
    end
  
    % Choose bisection or interpolation
    i_bis = find( (abs(e) < toler) | (abs(fa) <= abs(fb)) );
    i_int = find( ~( (abs(e) < toler) | (abs(fa) <= abs(fb)) ) );
    if ~isempty(i_bis)
        if trace > 3
            fprintf(fid,'#05# abs(e) < toler -> %d  -  abs(fa) <= abs(fb) -> %d\n', [(abs(e) < toler); (abs(fa) <= abs(fb))]);
        end
        % Bisection
        d(i_bis) = m(i_bis);  e(i_bis) = m(i_bis);
        procedure(i_bis) = 1;
        if trace > 3
            fprintf(fid,'#06# d = %13.6g  -  e = %13.6g\n', [d; e]);
        end
    end
    if ~isempty(i_int)
        % Interpolation
        procedure_int = procedure(i_int);
        e_int = e(i_int);
        d_int = d(i_int);
        s_int = fb(i_int)./fa(i_int);
        m_int = m(i_int);
        p_int = zeros(1,length(i_int));
        q_int = zeros(1,length(i_int));
        i_lin = find( a(i_int) == c(i_int) );
        i_quad = find( ~(a(i_int) == c(i_int)) );
        if trace > 3
            fprintf(fid,'#07# s_int = %13.6g\n', s_int);
        end
        if ~isempty(i_lin)
            % Linear interpolation
            p_int_lin = 2.0*m_int.*s_int;
            q_int_lin = 1.0 - s_int;
            p_int(i_lin) = p_int_lin(i_lin);
            q_int(i_lin) = q_int_lin(i_lin);
            if trace > 3
                fprintf(fid,'#08a# p_int = %13.6g  -  q_int = %13.6g\n', [p_int; q_int]);
            end
        end
        if ~isempty(i_quad)
            % Inverse quadratic interpolation
            q_int_quad = fa(i_int)./fc(i_int);
            r_int_quad = fb(i_int)./fc(i_int);
            p_int_quad = s_int.*(2.0*m_int.*q_int_quad.*(q_int_quad - r_int_quad) - (b(i_int) - a(i_int)).*(r_int_quad - 1.0));
            q_int_quad = (q_int_quad - 1.0).*(r_int_quad - 1.0).*(s_int - 1.0);
            p_int(i_quad) = p_int_quad(i_quad);
            q_int(i_quad) = q_int_quad(i_quad);
            if trace > 3
                fprintf(fid,'#08b# r_int_quad = %13.6g  -  p_int = %13.6g  -  q_int = %13.6g\n', [r_int_quad; p_int; q_int]);
            end
        end
        i_larger = find(p_int>0);
        i_lower = find(p_int<=0);
        if ~isempty(i_larger)
            q_int(i_larger) = -q_int(i_larger);
        end
        if ~isempty(i_lower)
            p_int(i_lower) = -p_int(i_lower);
        end
        if trace > 3
            fprintf(fid,'#09# p_int = %13.6g  -  q_int = %13.6g\n', [p_int; q_int]);
        end
        
        % Is interpolated point acceptable
        i_int_new = find( (2.0*p_int < 3.0*m_int.*q_int - abs(toler(i_int).*q_int)) & (p_int < abs(0.5*e(i_int).*q_int)) );
        i_bis_new = find( ~( (2.0*p_int < 3.0*m_int.*q_int - abs(toler(i_int).*q_int)) & (p_int < abs(0.5*e(i_int).*q_int)) ) );
        if ~isempty(i_int_new)
            e_int(i_int_new) = d_int(i_int_new); d_int(i_int_new) = p_int(i_int_new)./q_int(i_int_new);
            procedure_int(i_int_new) = 2;
            if trace > 3
                fprintf(fid,'#10a# d_int = %13.6g  -  e_int = %13.6g\n', [d_int; e_int]);
            end
        end
        if ~isempty(i_bis_new)
            d_int(i_bis_new) = m_int(i_bis_new);  e_int(i_bis_new) = m_int(i_bis_new);
            procedure_int(i_bis_new) = 1;
            if trace > 3
                fprintf(fid,'#10b# d_int = %13.6g  -  e_int = %13.6g\n', [d_int; e_int]);
            end
        end
        d(i_int) = d_int;  e(i_int) = e_int;
        procedure(i_int) = procedure_int;
    end
  
    % Next point
    a = b;
    fa = fb;
    if trace > 3
        fprintf(fid,'#11# a = %13.6g  -  fa = %13.6g\n', [a; fa]);
    end
    i_d = find( abs(d) > toler );
    if ~isempty(i_d)
        b(i_d) = b(i_d) + d(i_d);
        if trace > 3
            fprintf(fid,'#12a# b = %13.6g\n', b);
        end
    end
    i_bc = find( b > c & abs(d) <= toler );
    if ~isempty(i_bc)
        b(i_bc) = b(i_bc) - toler(i_bc);
        if trace > 3
            fprintf(fid,'#12b# b = %13.6g\n', b);
        end
    end
    i_else = find( abs(d) <= toler & b <= c );
    if ~isempty(i_else)
        b(i_else) = b(i_else) + toler(i_else);
        if trace > 3
            fprintf(fid,'#12c# b = %13.6g\n', b);
        end
    end

    fb = sigmafun(b,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
    fcount = fcount + 1*nx;

    iter = iter + 1;
    
    if trace > 3
        fprintf(fid,'#13# fb = %13.6g\n', fb);
        fprintf(fid,'#14# iter = %d\n', iter);
        fprintf(fid,'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\n');
    end

end % Main loop

% Output last chosen b
if trace > 2
    fprintf(fid,'Func-count = %5.0f\n', fcount);
    fprintf(fid,header2);
    fprintf(fid,'%13.6g %13.6g %13.6g\n',[b;fb;procedure]);
end

output.intervaliterations = intervaliter;
output.iterations = iter;
output.funcCount = fcount;
output.algorithm = 'bisection, interpolation';
fval = sigmafun(b,alla,rho,lsf,probdata,analysisopt,gfundata,femodel,randomfield);
fcount = fcount + 1*nx;

if all( abs(fval) <= max(abs(savefa),abs(savefb)) )
    msg = sprintf('Zero found in the interval [%g, %g]\n',[savea; saveb]);
else
    exitflag = -5; 
    msg = sprintf([...
        'Current point x may be near a singular point. The interval [%g, %g] \n', ... 
        'reduced to the requested tolerance and the function changes sign in the interval,\n', ...
        'but f(x) increased in magnitude as the interval reduced.\n'],[savea; saveb]);
end
if trace > 1
    fprintf(fid,'\n');
    fprintf(fid,msg);
end
output.message = msg;


if fid ~= 1
   fclose(fid);
end


%------------------------------------------------------------------

function [exitflag,msg] = disperr(y, fy, trace)
%DISPERR Display an appropriate error message when FY is Inf, 
%   NaN, or complex.  Assumes Y is the value and FY is the function 
%   value at Y. If FY is neither Inf, NaN, or complex, it generates 
%   an error message.

if ~isfinite(fy)  % NaN or Inf detected
    exitflag = -3;
    msg = ...
        sprintf(['Exiting fzero: aborting search for an interval containing a sign change\n' ...
                 '    because NaN or Inf function value encountered during search.\n' ...
                 '(Function value at %g is %g.)\n' ...
                 'Check function or try again with a different starting value.\n'],y,fy);
    if trace > 0
        fprintf(fid,msg);
    end
elseif ~isreal(fy) % Complex value detected
    exitflag = -4;
    msg = ...
        sprintf(['Exiting fzero: aborting search for an interval containing a sign change\n' ...
                 '    because complex function value encountered during search.\n' ...
                 '(Function value at %g is %s.)\n' ...
                 'Check function or try again with a different starting value.\n'],y,num2str(fy));
    if trace > 0
        fprintf(fid,msg);        
    end
else
    error('MATLAB:fzero:disperr:InvalidArg',...
        'DISPERR (in FZERO) called with invalid argument.')
end
