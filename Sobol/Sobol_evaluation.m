function [retour,retourT,retourTout] = Sobol_evaluation(lsf,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun

% Options settings

marg           = probdata.marg;
R              = probdata.correlation;
transf_type    = probdata.transf_type;
flag_sens      = probdata.flag_sens;
nrv = size(marg,1);

sampling       = analysisopt.sampling;
block_size     = analysisopt.block_size;
rand_generator = analysisopt.rand_generator;
num_sim        = analysisopt.num_sim;
first_indices  = analysisopt.first_indices;
total_indices  = analysisopt.total_indices;
all_indices    = analysisopt.all_indices;
Sobolsetopt    = analysisopt.Sobolsetopt;

if isfield(gfundata(lsf),'ng')
   ng = gfundata(lsf).ng;
else
   ng = 1;
end

% Initialization

k  = 0;  % Increment for block_size

S  = []; % First and all order indices
ST = []; % Total indices

f1  = zeros(1,ng);
s   = zeros(1,ng);
% f0  = 0; % Initialization for mean of y=f(x)=g(x)
% f00 = 0;
% D   = 0; % Initialization for variance of y=f(x)=g(x)
D1  = zeros(1,ng);
M   = [];% Matrix indexation of Sobol' indices: each ligne defines an indice

M      = create_matrix(nrv);
somme1 = zeros(size(M,1),ng);
somme2 = zeros(size(M,1),ng);
somme3 = zeros(size(M,1),ng);
summ   = zeros(size(M,1),ng);


% Simulation of samples

while k < num_sim

   block_size = min(block_size, num_sim-k );

   k = k + block_size;

   if sampling == 1  % Random sample for Sobol' indices assessment

      allu = []; % Matrix generated in normal space
      allx = []; % Matrix generated in physical space

      for r = 1:2
         switch rand_generator
            case 0
               allu (:,:,r) = randn(nrv,block_size);
            otherwise
% %                allu (:,:,r) = inv_norm_cdf(twister(nrv,block_size)); % original
               allu (:,:,r) = inv_norm_cdf(rand(nrv,block_size));
         end
         allx (:,:,r) = u_to_x(allu(:,:,r),probdata); % Transform into physical space
      end

      allx1 = allx (:,:,1);
      allx2 = allx (:,:,2);
      
      u = allx1';
      v = allx2';

   else % Quasi-random samples for Sobol' indices assessment
      
      allu = [];% matrix generated in normal space
      allv = [];% matrix generated in normal space
      allx = [];% matrix generated in original space
      al   = [];
      A    = [];

      if Sobolsetopt == 0 % Use of Statistical toolbox of Matlab with the sobolset function
         allp = sobolset(2);
         uv = inv_norm_cdf(allp((k-block_size+1+1):(k+1),:)');
      end

      if Sobolsetopt == 1 % Use of Sobol' sequence from BRODA
         al = zeros(2,block_size);
         for j = (k-block_size+1):k
            al(:,j) = sobolseq51(j,2)';
         end
         uv = inv_norm_cdf(al(:,(k-block_size+1):k));
      end

      % Quasi-random sequences by randomizing Sobol' sequences

      for i = 1:nrv
         permut    = randperm(block_size);
         utmp      = uv(1,:)';
         allu(:,i) = utmp(permut);
         permut    = randperm(block_size);
         utmp      = uv(2,:)';
         allv(:,i) = utmp(permut);
      end

      allx1 = u_to_x(allu',probdata);
      allx2 = u_to_x(allv',probdata);

      u = allx1';
      v = allx2';
      
   end

   numind = 0;       % Number of the Sobol' indices
   Di = []; DT = []; % Sobols'indices numerators initialization

   % g(u) and g(v)

   [ allg, dummy ] = gfun(lsf,u','no ',probdata,analysisopt,gfundata,femodel,randomfield); % Evaluate limit-state function
   yu = allg';
   
   [ allg, dummy ] = gfun(lsf,v','no ',probdata,analysisopt,gfundata,femodel,randomfield); % Evaluate limit-state function
   yv = allg';

   % First and second order moments assessment

   f1 = f1 + sum(yu,1);
   f0 = f1/k;          % Mean value

   D1 = D1 + sum(yu.^2,1);
   D  = D1/k - f0.^2;  % Variance

   % Modified estimation of the mean value by taking into account the two samples
   % and allowing a better estimation of Sobol' indices

   s   = s + sum(yu.*yv,1);
   f00 = s/k;

   % Assessment of Sobols' indices

   % Assessment of first indices
   
   for i = 2:size(M,1)
      
      liste = [];
      j = find(M(i,:)); % Allow to write the indexation number from matrix M
      ll = WriteIndex(j);
      numind = [ numind; ll ];

      if length(j) == 1 % For first Sobols' indices
         
         if first_indices == 1
            x_star = v;
            x_star(:,j) = u(:,j); % In order to calculate conditional moments (in "j")
            [ allg, dummy ] = gfun(lsf,x_star','no ',probdata,analysisopt,gfundata,femodel,randomfield);
            yv = allg'; % Second evaluation of g with v modified by "conditional to Xj"
            summ(i,:) = sum(yu.*yv,1);
            somme1(i,:) = somme1(i,:) + summ(i,:); % ûi*n
            Di(i,:) = somme1(i,:)/(k-1) - f00;
         end
         
      else

         % Assessment of all order (>1) indices

         if all_indices == 1 & nrv < 10
            x_star = v;
            for l = j
               x_star(:,l) = u(:,l); % "Conditional to several Xj"
            end
            [ allg, dummy ] = gfun(lsf,x_star','no ',probdata,analysisopt,gfundata,femodel,randomfield);
            yv = allg';
            summ(i,:) = sum(yu.*yv,1);
            somme2(i,:) = somme2(i,:) + summ(i,:);  % ûij*n
            Di(i,:) = somme2(i,:)/(k-1) - f00;
            liste = Inverse_create_matrix(j,nrv,M);
            for l = liste
               Di(i,:) = Di(i,:) - Di(l,:);
            end
         else
            Di(i,:) = 0;% Indices are not calculated
         end
         
      end

      if all_indices == 1 || first_indices == 1
         S(i,:) = Di(i,:)./D;
      end
      
   end

   % Assessment of total indices

   if total_indices == 1
      
      numindT = [];

      for i = 1:nrv
         x_star = u;
         x_star(:,i) = v(:,i);
         [ allg, dummy ] = gfun(lsf,x_star','no ',probdata,analysisopt,gfundata,femodel,randomfield);
         yxstar = allg';
         summ(i,:) = sum(yxstar.*yu,1);
         somme3(i,:) = somme3(i,:) + summ(i,:);
         DT(i,:) = somme3(i,:)/(k-1) - f00;
         numindT = [ numindT; i ];
         ST(i,:) = 1 - DT(i,:)./D;
      end
      
   end
   
end

% Results

retourT = []; retourTout = []; retour = []; % Results matrix initialization

if first_indices == 1 & all_indices == 0
   res        = [ numind  S ];
   res2       = sortrows(res,1);
   res        = res2(2:nrv+1,:);
   retour     = res;
end

if first_indices == 1 & all_indices == 1
   res        = [ numind  S ];   % Create table of Sobols' indices
   res2       = sortrows(res,1); % For ranking indices (first S1, second S2, S12, ...)
   res        = res2(2:end,:);
   retour     = res;
end

if total_indices == 1
   resT       = [ numindT  ST ];
   retourT    = sortrows(resT,1);
end

if all_indices == 1
   retourTout = [ retour; retourT ];
end


%--------------------------------------------------------------------
%              Functions used in Sobol_Evaluation.m
%--------------------------------------------------------------------

function M = create_matrix(nrv)
% Create a matrix from which each row defines an index
% Ex: the row 0101 represents the Sobol' index S_24
if nrv < 10
   N = 2^nrv;
   M = zeros(N,nrv);
   for k = 1:nrv
      v = [ zeros(2^(k-1),1); ones(2^(k-1),1) ];
      M(:,nrv-k+1) = repmat(v,2^(nrv-k),1);
   end
else
   % Identity matrix
   M = [ zeros(1,nrv); eye(nrv) ];
end

%--------------------------------------------------------------------

function x = fpart(xx)
x = xx - floor(xx);

%--------------------------------------------------------------------

function k = WriteIndex(j)
% Function for which 24 corresponds to [2 4]
k = 0;
for i = 1:length(j)
   k = k + 10^(length(j)-i)*j(i);
end

%--------------------------------------------------------------------

function liste = Inverse_create_matrix(indice,nrv,M)
% Inverse of CreateMatrix function,
% Ex: from 0101 you get [2 4]
liste = [];
bin = zeros(1,nrv);
for i = indice
   bin(i) = 1;
end
for i = 2:2^nrv
   d = bin - M(i,:);
   mauvais = 0;
   for j = 1:nrv
      if d(j) < 0, mauvais = 1; end
   end
   if mauvais == 0
      if sum(d)~=0, liste = [ liste i ]; end
   end
end