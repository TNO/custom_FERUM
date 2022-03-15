% Find misclassified points, if any

num_sim = length(allG(1,:));

Y = sign( allG(1,I_train) - g0 );
D = data(allU(:,I_train)',Y')

Gsvm = zeros(1,num_sim);

k = 0;
n_in_margin     = 0;
n_neg           = 0;
n_pos           = 0;
n_misclassified = 0;
block_size = analysisopt.buf_size;
svm_buf_size = analysisopt.svm_buf_size;
while k < num_sim

   block_size = min( block_size, num_sim-k );
   k = k + block_size;

   block_u = allU(:,(k-block_size+1):k);
   gridsel_option_svc = 0;
   block_g = eval_svm(block_u,svm_buf_size,SVC,gridsel_option_svc);

   Gsvm((k-block_size+1):k) = block_g;
   
   % Points in the margin
   I_in_margin = find( 1 > abs(block_g) );
   n_in_margin = n_in_margin + length(I_in_margin);

   I_neg = find( block_g <= 0 );
   n_neg = n_neg + length(I_neg);
   I_pos = find( block_g > 0 );
   n_pos = n_pos + length(I_pos);

   I_misclassified = find( ( allG(1,(k-block_size+1):k) - g0 ) .* block_g < 0 );
   n_misclassified = n_misclassified + length(I_misclassified);

end

fprintf(1,'\n');
fprintf(1,'Number of learning points            : %d\n', num_sim);
fprintf(1,'Number of negative points (real/SVM) : %d / %d\n', length( find( ( allG(1,:) - g0 ) <= 0 ) ), n_neg)
fprintf(1,'Number of positive points (real/SVM) : %d / %d\n', length( find( ( allG(1,:) - g0 ) > 0 ) ), n_pos)
fprintf(1,'Number of misclassified points       : %d\n', n_misclassified)
fprintf(1,'\n');

svmmisclassified(1,kk) = num_sim;
svmmisclassified(2,kk) = n_neg;
svmmisclassified(3,kk) = n_pos;
svmmisclassified(4,kk) = n_misclassified;
