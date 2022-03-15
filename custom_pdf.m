% Probability distrtibution function of a custom random variable
%
% p = CUSTOM_PDF(x, ID, type)
%
%INPUT
% x      - value at which we are interested in the pdf value,
% ID     - ID value of the distribution (specified in probdata.marg)
% type   - definition type of the custom distribution, 'sample' or 'point'
%
%OUTPUT
% p      - probability density value, pdf(X=x)
%
%NOTE:
%

function p = custom_pdf(x, ID, type)

persistent S

switch (lower(type))
    
    case 'sample'
        load(['tmp\kernel_sample_', num2str(ID), '.mat'], 'x_grid', 'pdf')
        p = interp1(x_grid, pdf, x);
        %p = interp1(x_grid, pdf, x, 'pchip');
    case 'point'
        % distribution function is given by discrete points
        if isempty(S)
            load('tmp\vector_distr_struct.mat', 'S')
        end
        
        x_grid = S(ID).x_grid;
        pdf    = S(ID).pdf;
        
        p = interp1(x_grid, pdf, x);
        
% %         % an attempt to avoid numerical problems
% %         tpdf        = norminv(pdf);
% %         idx         = ~isinf(tpdf);
% %         tpdf        = tpdf(idx);
% %         tx_grid     = x_grid(idx);
% % 
% %         p           = normcdf(interp1(tx_grid, tpdf, x));
end

end





