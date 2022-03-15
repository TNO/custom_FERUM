% Define here the directory where FERUM is installed
s = 'D:\Home\work\FERUM4.1';
addpath(genpath(s));

% twister('state',100*clock);

inputfile_Dakessian_example1

analysisopt.echo_flag            = 0;
analysisopt.analysistype         = 23;
analysisopt.flag_cov_pf_bounds   = 1;
analysisopt.block_size           = 500000;

num_sim = [ 1000 10000 ];

ntimes = 500;

pf              = zeros(length(num_sim),ntimes);
cov_pf          = zeros(length(num_sim),ntimes);
cov_pf_bornemax = zeros(length(num_sim),ntimes);
beta            = zeros(length(num_sim),ntimes);
nfun            = zeros(length(num_sim),ntimes);

for j = 1:length(num_sim)
        
   analysisopt.num_sim = num_sim(j);
   
   k = 0;
   while k < ntimes
      k = k + 1;
      fprintf('j = %d - k = %d\n',j,k);
      ferum
      pf(j,k) = subsetsimulationresults.pf;
      cov_pf(j,k) = subsetsimulationresults.cov_pf;
      cov_pf_bornemax(j,k) = subsetsimulationresults.SubsetData.cov_pf_bounds(2);
      beta(j,k) = subsetsimulationresults.beta;
      nfun(j,k) = subsetsimulationresults.nfun;
      clear subsetsimulationresults
   end

end

for j = 1:length(num_sim)
   PF_MEAN(j) = mean(pf(j,:));
   COV_PF_EMP(j) = std(pf(j,:)) / mean(pf(j,:));
   COV_PF(j) = mean(cov_pf(j,:));
   COV_PF_BORNEMAX(j) = mean(cov_pf_bornemax(j,:));
   NFUN(j) = mean(nfun(j,:));
end
PF_MEAN
COV_PF_EMP
COV_PF
COV_PF_BORNEMAX