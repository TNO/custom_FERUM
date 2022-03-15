% Convert bounds contrains to linear constrains for NLopt COBYLA

function opt = cobyla_bound_con(lb, ub, opt)

if nargin < 3
    opt = [];
end

idx_lb  = find(lb ~= -Inf);
n_lb    = length(idx_lb);

idx_ub  = find(ub ~= Inf);
n_ub    = length(idx_ub);

n_tot   = n_lb + n_ub;
fc      = cell(1, n_tot);

for ii = 1:n_lb
    loc     = ii;
    x_ii    = idx_lb(ii);
    fc{loc} = eval(['@(x) ', sprintf('%.16f', lb(x_ii)), ' - x(', num2str(x_ii), ')']);
end
for ii = 1:n_ub
    loc     = n_lb + ii;
    x_ii    = idx_ub(ii);
    fc{loc} = eval(['@(x) ', ' x(',num2str(x_ii),') - ', sprintf('%.16f', ub(x_ii)) ]);
end

opt.fc        = fc;
opt.fc_tol    = 1e-4*ones(1, n_tot);

end