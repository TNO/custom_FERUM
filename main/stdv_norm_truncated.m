function stdv_t = stdv_norm_truncated(mean,stdv,xmin,xmax)

mean_t = mean_norm_truncated(mean,stdv,xmin,xmax);

var_t = 1/(normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) * ...
           ( mean^2 * (normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) - ...
             2*mean*stdv * (normpdf((xmax-mean)/stdv)-normpdf((xmin-mean)/stdv)) - ...
             stdv^2 * ((xmax-mean)/stdv*normpdf((xmax-mean)/stdv)-(xmin-mean)/stdv*normpdf((xmin-mean)/stdv)) + ...
             stdv^2 * (normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) ) - ...
        mean_t^2;

stdv_t = var_t^0.5;