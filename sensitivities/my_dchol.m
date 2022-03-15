  function [L,dL,ierr] = my_dchol(A,dA);
   %   
      [n,n] = size(A);
      ierr = 0;
   %

      dLsurdA = dA; 
    
      for k = 1:n, 
   %
   %     exit if A is not positive definite
   %
         if (A(k,k) <= 0), ierr = k; fprintf('Error in Choleski decomposition - Matrix must be positive definite\n'); return; end
   %
   %     Compute main diagonal elt and then scale the k-th column 
   %
         dLsurdA(k,k) = 1/2 * dLsurdA(k,k) * 1/sqrt(A(k,k)) ;
         A(k,k) = sqrt(A(k,k));

         dLsurdA(k+1:n,k) = ( dLsurdA(k+1:n,k) * A(k,k) - A(k+1:n,k) * dLsurdA(k,k) ) / A(k,k)^2 ;
         A(k+1:n,k) = A(k+1:n,k)/A(k,k);
   %
   %     Update lower triangle of the trailing (n-k) by (n-k) block
   %
         for j = k+1:n,
            dLsurdA(j:n,j) = dLsurdA(j:n,j) - dLsurdA(j:n,k) * A(j,k) - A(j:n,k) * dLsurdA(j,k);
            A(j:n,j) = A(j:n,j) - A(j:n,k)*A(j,k);
         end
         
      end

      L = tril(A);
      dL = tril(dLsurdA);