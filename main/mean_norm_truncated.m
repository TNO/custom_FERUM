function mean_t = mean_norm_truncated(mean,stdv,xmin,xmax)

mean_t = 1/(normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) * ...
            ( mean * (normcdf((xmax-mean)/stdv)-normcdf((xmin-mean)/stdv)) - ...
              stdv * (normpdf((xmax-mean)/stdv)-normpdf((xmin-mean)/stdv)) );