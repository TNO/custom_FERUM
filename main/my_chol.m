  function [L,ierr] = my_chol(A);
   %   
      [n,n] = size(A);
      ierr = 0;
   %

      for k = 1:n
   %
   %     exit if A is not positive definite
   %
         if (A(k,k) <= 0), ierr = k; fprintf('Error in Choleski decomposition - Matrix must be positive definite\n'); return; end
   %
   %     Compute main diagonal elt. and then scale the k-th column 
   %
         A(k,k) = sqrt(A(k,k));
         A(k+1:n,k) = A(k+1:n,k)/A(k,k);
   %
   %     Update lower triangle of the trailing (n-k) by (n-k) block
   %
         for j = k+1:n,
            A(j:n,j) = A(j:n,j) - A(j:n,k)*A(j,k);
         end

      end

      L = tril(A);