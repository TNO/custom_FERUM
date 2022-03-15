function auxva_string = auxva(VA,expx)

% Called by gen_replace_awk.m.
% Show how to insert a string (which may be composed of multiple lines), result of a function
% of the current realizations of the random variables.
% The field femodel.auxva must exist and set to 1, if we want gen_replace_awk.m to call auxva_string.m.
%
% Nota: This file is not specific to SSLL116b example.

for iexpx = 1:size(expx,1)    
   x = expx(iexpx,1);
   eval([ VA(iexpx).name ' = expx(' num2str(iexpx) ',1);']);
end


sigma = 25:25:1250;
eps = sigma/E + alpha0*sigmay/E * (sigma/sigmay).^n0;
auxva_string1 = '';
for i = 1:(length(sigma)-1)
   auxva_string1 = [ auxva_string1 sprintf(' %.10e , %.10e ,\\n', [eps(i); sigma(i)]) ];
end
auxva_string1 = [ auxva_string1 sprintf(' %.10e , %.10e ,', [eps(length(sigma)); sigma(length(sigma))]) ];

auxva(2) = P * Ri^2 / ( (Ri+t)^2 - Ri^2 );
auxva_string2 = sprintf('%.10e', auxva(2));

auxva(3) = sigmat + P * Ri^2 / ( (Ri+t)^2 - Ri^2 );
auxva_string3 = sprintf('%.10e', auxva(3));

auxva_string = { auxva_string1 auxva_string2 auxva_string3 };
