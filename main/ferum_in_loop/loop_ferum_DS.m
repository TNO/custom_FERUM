% Define here the directory where FERUM is installed
s = 'D:\Home\work\FERUM4.1';
addpath(genpath(s));

% twister('state',100*clock);

inputfile_Dakessian_example1

analysisopt.echo_flag            = 0;
analysisopt.analysistype         = 22;
analysisopt.block_size           = 500000;

analysisopt.dir_flag             = 'random'; % 'det': deterministic points uniformly distributed on the unit hypersphere using eq_point_set.m function
                                             % 'random': random points uniformly distributed on the unit hypersphere

num_sim = [ 50 100 500 1000 5000 10000 ];

ntimes = 1000;

pf     = zeros(length(num_sim),ntimes);
cov_pf = zeros(length(num_sim),ntimes);
beta   = zeros(length(num_sim),ntimes);
nfun   = zeros(length(num_sim),ntimes);

for j = 1:length(num_sim)
        
   analysisopt.num_sim = num_sim(j);
   
   k = 0;
   while k < ntimes
      k = k + 1;
      fprintf('j = %d - k = %d\n',j,k);
      ferum
      pf(j,k) = dirsimulationresults.pf;
      cov_pf(j,k) = dirsimulationresults.cov_pf;
      beta(j,k) = dirsimulationresults.beta;
      nfun(j,k) = dirsimulationresults.nfun;
      clear dirsimulationresults
   end

end

for j = 1:length(num_sim)
   PF_MEAN(j) = mean(pf(j,:));
   COV_PF_EMP(j) = std(pf(j,:)) / mean(pf(j,:));
   COV_PF(j) = mean(cov_pf(j,:));
   NFUN(j) = mean(nfun(j,:));
end
PF_MEAN
COV_PF_EMP
COV_PF