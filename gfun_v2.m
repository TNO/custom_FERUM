function [ G, grad_g ] = gfun_v2(lsf,x,grad_flag,probdata,analysisopt,gfundata,femodel,randomfield)

global nfun

% Takes a default value for analysisopt.ffdpara, if not defined in input file
if ~isfield(analysisopt,'ffdpara')
    switch lower(gfundata(lsf).evaluator)
        case 'basic'
            analysisopt.ffdpara = 1000;
        otherwise
            analysisopt.ffdpara = 50;
    end
end


nrv = size(x,1);
nx  = size(x,2);
if isfield(gfundata(lsf),'thetag')
    nthetag = size(gfundata(lsf).thetag,2);
else
    nthetag = 0;
end


% If multi_proc option not specified for 'basic' limit-state function, the g-function
% is supposed to be of non-vectorized type and a sequential process is forced
if ~isfield(analysisopt,'multi_proc') && gfundata(lsf).evaluator == 'basic'
    analysisopt.multi_proc = 0;
end


switch analysisopt.multi_proc
     
    % sequential calls
    case 0
        
        if strcmp(grad_flag, 'no')
            
            if nx > 1
                
                if isfield(gfundata(lsf),'ng')
                    ng = gfundata(lsf).ng;
                else
                    ng = 1;
                end
                G = zeros(ng,nx);
                grad_g = zeros(1,nx);
                for i = 1:nx
                    switch lower(gfundata(lsf).evaluator)
                        case 'basic'
                            [ g, dummy ] = gfunbasic_v2(lsf,x(:,i),'yes','no',probdata,analysisopt,gfundata);
                            nfun = nfun+1;
                        otherwise
                            eval(['[ g, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x(:,i),''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                    end
                    if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                        g = gfunwithbulge(lsf,x(:,i),probdata,gfundata,g);
                    end
                    G(:,i) = g;
                    grad_g(i) = dummy;
                end
                
            elseif nthetag > 1
                
                G = zeros(1,nthetag);
                original_thetag = gfundata(lsf).thetag;
                for i = 1:nthetag
                    gfundata(lsf).thetag = original_thetag(:,i);
                    switch lower(gfundata(lsf).evaluator)
                        case 'basic'
                            [ g, ~ ] = gfunbasic_v2(lsf,x,'yes','no',probdata,analysisopt,gfundata);
                            nfun = nfun+1;
                        otherwise
                            eval(['[ g, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                    end
                    G(i) = g;
                end
                grad_g = 0;
                
            else
                
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        [ G, dummy ] = gfunbasic_v2(lsf,x,'yes','no',probdata,analysisopt,gfundata);
                        nfun = nfun+1;
                    otherwise
                        eval(['[ G, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                end
                if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                    G = gfunwithbulge(lsf,x,probdata,gfundata,G);
                end
                grad_g = dummy;
                
            end
        % <---> 
        elseif strcmp(grad_flag, 'ffd')
            
            allx = x;
            G = zeros(1,nx);
            grad_g = zeros(nrv,nx);
            
            for i = 1:nx
                
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        [ g, ~ ] = gfunbasic_v2(lsf,allx(:,i),'yes','no',probdata,analysisopt,gfundata);
                        nfun = nfun+1;
                    otherwise
                        eval(['[ g, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,allx(:,i),''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                end
                if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                    g = gfunwithbulge(lsf,allx(:,i),probdata,gfundata,g);
                end
                G(i) = g;
                
                marg = probdata.marg;
                original_x = allx(:,i);
                for j = 1:nrv
                    x = original_x;
                    h = marg(j,3)/analysisopt.ffdpara;
                    x(j) = x(j) + h;
                    switch lower(gfundata(lsf).evaluator)
                        case 'basic'
                            [ g_a_step_ahead, ~ ] = gfunbasic_v2(lsf,x,'yes','no',probdata,analysisopt,gfundata);
                            nfun = nfun+1;
                        otherwise
                            eval(['[ g_a_step_ahead, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                    end
                    if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                        g_a_step_ahead = gfunwithbulge(lsf,x,probdata,gfundata,g_a_step_ahead);
                    end
                    grad_g(j,i) = (g_a_step_ahead - g)/h;
                end
                
            end
            
        elseif strcmp(grad_flag, 'ddm')
            
            switch lower(gfundata(lsf).evaluator)
                case 'basic'
                    [ G, grad_g ] = gfunbasic_v2(lsf,x,'yes','yes',probdata,analysisopt,gfundata);
                    nfun = nfun+1;
                otherwise
                    eval(['[ G, grad_g ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''yes'',probdata,analysisopt,gfundata,femodel,randomfield);']);
            end
            
        % <--->    
        % =====================================================================================================================================================
        % EXTENSION - Arpad Rozsas, 2017-April-18
        % mixed 'ffd' and 'ddm' calculation of the gradient
        %
        % add check of input! if 'ffd_ddm' is selected!!
        % what gfunwithbulge is doing, left intact for now
        elseif strcmp(grad_flag, 'ffd_ddm')

            allx = x;
            G       = nan(1,nx);
            grad_g  = nan(nrv,nx);
            
            for i = 1:nx
                
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        [ g, ~ ] = gfunbasic_v2(lsf,allx(:,i),'yes','no',probdata,analysisopt,gfundata);
                        nfun = nfun+1;
                    otherwise
                        eval(['[ g, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,allx(:,i),''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                end
                if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                    g = gfunwithbulge(lsf,allx(:,i),probdata,gfundata,g);
                end
                G(i) = g;
                
                marg = probdata.marg;
                original_x = allx(:,i);
                gradient_expression = gfundata(lsf).dgdq;
                for j = 1:nrv
                    x = original_x;
                    
                    % use 'ffd' to estimate the partial derivative
                    if any(isnan(gradient_expression{j}))
                        h = marg(j,3)/analysisopt.ffdpara;
                        x(j) = x(j) + h;
                        switch lower(gfundata(lsf).evaluator)
                            case 'basic'
                                [ g_a_step_ahead, ~ ] = gfunbasic_v2(lsf,x,'yes','no',probdata,analysisopt,gfundata);
                                nfun = nfun+1;
                            otherwise
                                eval(['[ g_a_step_ahead, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                        end
                        if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                            g_a_step_ahead = gfunwithbulge(lsf,x,probdata,gfundata,g_a_step_ahead);
                        end
                        grad_g(j,i) = (g_a_step_ahead - g)/h;
                        
                    % use 'ddm' with the user specified partial derivative
                    % only the gradient is needed as G is already
                    % available
                    else
                        % evaluates the gradient multiple times
                        % not elegant, but it is cheap with the current
                        % formulation, later might be extended, if the
                        % gradient is expensive
                        [ ~, grad_g_tmp ] = gfunbasic_v2(lsf,x,'no','yes',probdata,analysisopt,gfundata);
                        
                        logi_idx_exp  = ~cellfun(@(C) isnumeric(C) && any(isnan(C(:))), gfundata(lsf).dgdq);
                        grad_g(logi_idx_exp) = grad_g_tmp(logi_idx_exp);
                    end
                end
                
            end
        % =====================================================================================================================================================
        % <--->        
        else
            
            disp('ERROR: Invalid method for gradient computations');
            
        end
        
    % simultaneous calls    
    case 1
        
        block_size = analysisopt.block_size;
        
        if strcmp(grad_flag,'no')
            
            if nx > 1
                
                if isfield(gfundata(lsf),'ng')
                    ng = gfundata(lsf).ng;
                else
                    ng = 1;
                end
                G = zeros(ng,nx);
                dummy = zeros(1,nx);
                
                k = 0;
                while k < nx
                    
                    block_size = min( block_size, nx-k );
                    blockx = x(:,(k+1):(k+block_size));
                    
                    switch lower(gfundata(lsf).evaluator)
                        case 'basic'
                            [ blockG, blockdummy ] = gfunbasic_v2(lsf,blockx,'yes','no',probdata,analysisopt,gfundata);
                        otherwise
                            eval(['[ blockG, blockdummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,blockx,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                    end
                    if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                        blockG = gfunwithbulge(lsf,blockx,probdata,gfundata,blockG);
                    end
                    
                    G(:,(k+1):(k+block_size)) = blockG;
                    dummy(1,(k+1):(k+block_size)) = blockdummy;
                    
                    k = k + block_size;
                    
                end
                
                grad_g = dummy;
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        nfun = nfun+nx;
                end
                
            elseif nthetag > 1
                
                G = zeros(1,nthetag);
                dummy = zeros(1,nthetag);
                
                k = 0;
                while k < nthetag
                    
                    block_size = min( block_size, nthetag-k );
                    blockgfundata = gfundata;
                    blockgfundata(lsf).thetag = blockgfundata(lsf).thetag(:,(k+1):(k+block_size));
                    
                    switch lower(gfundata(lsf).evaluator)
                        case 'basic'
                            [ blockG, blockdummy ] = gfunbasic_v2(lsf,x,'yes','no',probdata,analysisopt,blockgfundata);
                        otherwise
                            eval(['[ blockG, blockdummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''no'',probdata,analysisopt,blockgfundata,femodel,randomfield);']);
                    end
                    
                    G(1,(k+1):(k+block_size)) = blockG;
                    dummy(1,(k+1):(k+block_size)) = blockdummy;
                    
                    k = k + block_size;
                    
                end
                
                grad_g = dummy;
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        nfun = nfun+nthetag;
                end
                
            else
                
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        [ G, dummy ] = gfunbasic_v2(lsf,x,'yes','no',probdata,analysisopt,gfundata);
                    otherwise
                        eval(['[ G, dummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                end
                if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                    G = gfunwithbulge(lsf,x,probdata,gfundata,G);
                end
                
                grad_g = dummy;
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        nfun = nfun+1;
                end
                
            end
            
        elseif strcmp(grad_flag, 'ffd')
            
            allx = zeros(nrv,nx*(1+nrv));
            allx(:,1:(1+nrv):(1+(nx-1)*(1+nrv))) = x;
            allh = zeros(1,nrv);
            
            marg = probdata.marg;
            original_x = x;
            for j = 1:nrv
                x = original_x;
                allh(j) = marg(j,3)/analysisopt.ffdpara;
                x(j,:) = x(j,:) + allh(j)*ones(1,nx);
                allx(:,(1+j):(1+nrv):(1+j+(nx-1)*(1+nrv))) = x;
            end
            
            allG = zeros(1,nx*(1+nrv));
            
            k = 0;
            while k < nx*(1+nrv)
                
                block_size = min( block_size, nx*(1+nrv)-k );
                blockx = allx(:,(k+1):(k+block_size));
                
                switch lower(gfundata(lsf).evaluator)
                    case 'basic'
                        [ blockG, blockdummy ] = gfunbasic_v2(lsf,blockx,'yes','no',probdata,analysisopt,gfundata);
                    otherwise
                        eval(['[ blockG, blockdummy ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,blockx,''no'',probdata,analysisopt,gfundata,femodel,randomfield);']);
                end
                if isfield(gfundata(lsf),'bulge') && ( gfundata(lsf).bulge == 1 )
                    blockG = gfunwithbulge(lsf,blockx,probdata,gfundata,blockG);
                end
                
                allG(1,(k+1):(k+block_size)) = blockG;
                
                k = k + block_size;
                
            end
            
            G = allG(1:(1+nrv):(1+(nx-1)*(1+nrv)));
            grad_g = zeros(nrv,nx);
            for j = 1:nrv
                grad_g(j,:) = ( allG((1+j):(1+nrv):((1+j)+(nx-1)*(1+nrv))) - G ) / allh(j);
            end
            
            switch lower(gfundata(lsf).evaluator)
                case 'basic'
                    nfun = nfun+nx*(1+nrv);
            end
            
        elseif grad_flag == 'ddm'
            
            switch lower(gfundata(lsf).evaluator)
                case 'basic'
                    [ G, grad_g ] = gfunbasic(lsf,x,'yes',probdata,analysisopt,gfundata);
                otherwise
                    eval(['[ G, grad_g ] = gfun' lower(gfundata(lsf).evaluator) '(lsf,x,''yes'',probdata,analysisopt,gfundata,femodel,randomfield);']);
            end
            
            switch lower(gfundata(lsf).evaluator)
                case 'basic'
                    nfun = nfun+1;
            end
            
        end
        
    otherwise
        
        disp('multi_proc option of analysisopt incorrectly defined!');
        
        
end
