function block_g = eval_svm(block_u,svm_buf_size,Separator,gridsel_option)

num_sim = size(block_u,2);

yEstim_data = data(block_u',zeros(1,num_sim)');

if gridsel_option == 1
   sz = 1;
else
   sz = get_dim(Separator.Xsv);
   if sz == 0, sz = 1; end
end

% fprintf('num_sim = %d - sz = %d - dim(yEstim_data) = %d\n', num_sim, sz, get_dim(yEstim_data));

batch_size = round(svm_buf_size/sz);
if get_dim(yEstim_data) > batch_size
%    fprintf(1,'Nb of loops required for eval_svm: %d\n', round(get_dim(yEstim_data)/batch_size));
end

block_g = zeros(1,get_dim(yEstim_data));

for i = 1:batch_size:get_dim(yEstim_data)
   take = i:min(i+batch_size-1,get_dim(yEstim_data));
   if ~isempty(Separator.alpha)
      kerMaTemp = get_kernel(Separator.child,get(yEstim_data,take),Separator.Xsv);
      block_g(take) = Separator.alpha'* kerMaTemp + Separator.b0;
   else
      block_g(take) = Separator.b0*ones(1,length(take));
   end
end