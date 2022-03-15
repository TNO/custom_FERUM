function data = subset_y_threshold(data,analysisopt)

% Evaluate the threshold value y of the intermediate limit-state-function
% Make preliminary calculations for estimating the coefficient of variation of the failure probability pf

p = data.pf_target;
  
IndInf = data.Indices(end,1);
IndSup = data.Indices(end,2);

G = data.G(IndInf:IndSup);
N = length(G);

y = ferum_quantile(G,p);

if y < 0
   y = 0;
   p = length( find( G < y ) ) / N;
end

data.y = [ data.y y ];
data.p = [ data.p p ];

Indgerm = find( G < y );
Indgerm = Indgerm + IndInf - 1;

if data.Nb_step == 0

   data.Indgerm = Indgerm;

else

   if size(Indgerm,2) == size(data.Indgerm,2)
      data.Indgerm = [ data.Indgerm; Indgerm ];
   elseif size(Indgerm,2) > size(data.Indgerm,2)
      data.Indgerm = [ data.Indgerm zeros(size(data.Indgerm,1),size(Indgerm,2)-size(data.Indgerm,2)) ];
      data.Indgerm = [ data.Indgerm; Indgerm ];
   elseif size(Indgerm,2) < size(data.Indgerm,2),
      data.Indgerm = [ data.Indgerm; [Indgerm zeros(1,size(data.Indgerm,2)-size(Indgerm,2))] ];
   end

end


if isfield(analysisopt,'flag_cov_pf_bounds') & analysisopt.flag_cov_pf_bounds == 1
   
   if data.Nb_step == 0

      cov_pf_step = sqrt( (1-p) / (p*N) );

   else

      Nc = length(find(data.Indgerm(end-1,:)));
      N = IndSup - IndInf + 1;
      
% tic
% 
%       Itemp = zeros(size(G'));
%       Itemp( find( G <= y ) ) = 1;
%       I = [];
%       for i = 1:N/Nc
%          I = [ I Itemp((i-1)*Nc+1:i*Nc) ];
%       end
% 
%       Rzero = p*(1-p);
%       R = [];
%       for k = 1:(N/Nc-1)
%          R(k) = 0;
%          for j = 1:Nc
%             ind = find( I(j,:) );
%             indou = find( ind <= (N/Nc-k) );
%             R(k) = R(k) + sum( I(j,ind(indou)) .* I(j,ind(indou)+k) );
%          end
%          R(k) = 1/(N-k*Nc) * R(k) - p^2;
%       end
% 
%       rho = []; % rho(k) = R(k)/Rzero
%       gamma = 0;
%       for k = 1:(N/Nc-1)
%          rho(k) = R(k) / Rzero;
%          gamma = gamma + 2 * (1-k*Nc/N) * rho(k);
%       end
% 
%       cov_pf_step = sqrt( (1-p) / (p*N) * (1+gamma) );
% 
% toc

      Rzero = p*(1-p);
      
      I =  [ G <= y ];
      Indices = reshape((1:N),[],N/Nc)';
      gamma = 0;
      for k = 1:(N/Nc-1)
         Z = 0;
         for j = 1:Nc
            for l = 1:(N/Nc-k)
               Z = Z + I(Indices(l,j)) * I(Indices(l+k,j));
            end
         end
         rho(k) = ( 1/(N-k*Nc) * Z - p^2 ) / Rzero;
         gamma = gamma + 2 * (1-k*Nc/N) * rho(k);
      end

      cov_pf_step = sqrt( (1-p) / (p*N) * (1+gamma) );

   end
   
   data.cov_pf_step = [ data.cov_pf_step cov_pf_step ];
   
end