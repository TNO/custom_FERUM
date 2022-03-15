function formresults = postprocess_cobyla(dname, formresults)

fid     = fopen(dname, 'r');
dia       = fscanf(fid, '%c');

% ........................................................
% Beta values
% ........................................................
kyw     = 'nlopt_optimize eval';
idx     = strfind(dia, kyw);
idx1    = idx;
idx2    = idx + length(kyw) + 20;

nfun    = length(idx);

Recorded_beta = nan(nfun+1, 1);
for ii = 1:nfun
    
    tmp = sscanf(dia(idx1(ii):idx2(ii)), 'nlopt_optimize eval #%i: %f');
    
    Recorded_beta(ii) = tmp(2);
end

Recorded_beta(end) = formresults.beta;

% ........................................................
% Running times
% ........................................................
kyw     = 'Job has completed!';
idx     = strfind(dia, kyw);
idx1    = idx;
idx2    = idx + length(kyw) + 20;

n       = length(idx);

gfun_eval_time = nan(n, 1);
for ii = 1:n    
    tmp = sscanf(dia(idx1(ii):idx2(ii)), 'Job has completed! (%f sec)');
    gfun_eval_time(ii) = tmp(1);
end

% ........................................................
% gfun values
% ........................................................
kyw     = 'g        =';
idx     = strfind(dia, kyw);
idx1    = idx;
idx2    = idx + length(kyw) + 20;

n       = length(idx);

Recorded_g = nan(n, 1);
for ii = 1:n    
    tmp = sscanf(dia(idx1(ii):idx2(ii)), 'g        = %f');
    Recorded_g(ii) = tmp(1);
end

% ---------------------------------------------------------
% GATHER results
% ---------------------------------------------------------
% the last is not printed to the command window
formresults.nfun                    = nfun+1; 
formresults.Recorded_beta           = Recorded_beta;
formresults.Recorded_g              = Recorded_g;
formresults.gfun_eval_time          = gfun_eval_time;
formresults.gfun_eval_total_time    = sum(gfun_eval_time);
formresults.diary                   = dia;

end