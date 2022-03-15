function [probdata,gfundata,analysisopt] = update_data(lsf,probdata,analysisopt,gfundata,femodel)

% clear potential persistent variables
clear custom_pdf custom_cdf custom_invcdf
clear nonparametric_pdf nonparametric_cdf nonparametric_invcdf

% Assigns names to probdata.marg variables if not defined
marg = probdata.marg;
if ~isfield(probdata, 'name')
    eval([ 'margname = {' sprintf(' ''x%d''',[1:size(marg,1)]) ' };']);
else
    margname = probdata.name;
end
if size(margname,1) > size(margname,2)
    margname = margname';
end
correlation = probdata.correlation;


% Assigns names to gfundata(lsf).thetag variables if not defined
if isfield(gfundata(lsf), 'thetag')
    thetag = gfundata(lsf).thetag;
    if size(thetag,1) < size(thetag,2)
        thetag = thetag';
    end
    if ~isfield(gfundata(lsf), 'thetagname')
        eval([ 'thetagname = {' sprintf(' ''thetag%d''',[1:length(thetag)]) ' };']);
    else
        thetagname = gfundata(lsf).thetagname;
    end
    if size(thetagname,1) > size(thetagname,2)
        thetagname = thetagname';
    end
else
    thetagname = [];
    thetag = [];
end


% Assigns names to gfundata(lsf).cg variables if not defined
if isfield(gfundata(lsf), 'cg')
    cg = gfundata(lsf).cg;
    if size(cg,1) < size(cg,2)
        cg = cg';
    end
    if ~isfield(gfundata(lsf), 'cgname')
        eval([ 'cgname = {' sprintf(' ''cg%d''',[1:length(cg)]) ' };']);
    else
        cgname = gfundata(lsf).cgname;
    end
    if size(cgname,1) > size(cgname,2)
        cgname = cgname';
    end
else
    cgname = [];
    cg = [];
end


% Looks for redundant parameters in gfundata(lsf).thetag and probdata.marg
% Eliminates those in gfundata(lsf).thetag
if ~isempty(thetag)
    Itoremove = [];
    for i = 1:length(thetag)
        if ~isempty(strmatch(thetagname(i),margname,'exact'))
            Itoremove = [ Itoremove i ];
        end
    end
    if ~isempty(Itoremove)
        thetag(Itoremove)     = [];
        thetagname(Itoremove) = [];
    end
end


% Looks for redundant parameters in gfundata(lsf).cg and probdata.marg
% Eliminates those in gfundata(lsf).cg
if ~isempty(cg)
    Itoremove = [];
    for i = 1:length(cg)
        if ~isempty(strmatch(cgname(i),margname,'exact'))
            Itoremove = [ Itoremove i ];
        end
    end
    if ~isempty(Itoremove)
        cg(Itoremove)     = [];
        cgname(Itoremove) = [];
    end
end


switch lower(gfundata(lsf).evaluator)
    
    case 'basic'
        
        % Identifies r.v.'s, reliability parameters and deterministic parameters in probdata.marg
        Ix      = find(marg(:,1)>0);
        Ithetag = find(marg(:,1)==-1);
        Icg     = find(marg(:,1)==0);
        
        % Looks for reliability parameters in probdata.marg and transfers them to gfundata(lsf).thetag
        if ~isempty(Ithetag)
            thetagname = [ thetagname margname(Ithetag) ];
            thetag     = [ thetag; marg(Ithetag,2) ];
        end
        
        % Looks for deterministic parameters in probdata.marg and transfers them to gfundata(lsf).cg
        if ~isempty(Icg)
            cgname = [ cgname margname(Icg) ];
            cg     = [ cg; marg(Icg,2) ];
        end
        
        margname    = margname(Ix);
        marg        = marg(Ix,:);
        correlation = correlation(Ix,Ix);
        
        if isfield(gfundata(lsf), 'dgdq')
            gfundata(lsf).dgdq = gfundata(lsf).dgdq(Ix);
        end
        
    otherwise
        
        if ~isempty(femodel)
            
            data     = femodel.data;
            dataname = femodel.dataname;
            if size(dataname,1) > size(dataname,2)
                dataname = dataname';
            end
            
            [data,I] = sortrows(data,2);
            dataname = dataname(I);
            
            % Looks for redundant parameters in gfundata(lsf).thetag and femodel.data
            % Eliminates those in gfundata(lsf).thetag
            if ~isempty(thetag)
                Itoremove = [];
                for i = 1:length(thetag)
                    if ~isempty(strmatch(thetagname(i),dataname,'exact'))
                        Itoremove = [ Itoremove i ];
                    end
                end
                if ~isempty(Itoremove)
                    thetag(Itoremove)     = [];
                    thetagname(Itoremove) = [];
                end
            end
            
            % Looks for redundant parameters in gfundata(lsf).cg and femodel.data
            % Eliminates those in gfundata(lsf).cg
            if ~isempty(cg)
                Itoremove = [];
                for i = 1:length(cg)
                    if ~isempty(strmatch(cgname(i),dataname,'exact'))
                        Itoremove = [ Itoremove i ];
                    end
                end
                if ~isempty(Itoremove)
                    cg(Itoremove)     = [];
                    cgname(Itoremove) = [];
                end
            end
            
            % Updates probdata.marg and femodel.data
            for i = 1:size(data,1)
                % Forces type of variable in probdata.marg according to femodel.data
                if ~isnan(data(i,2))
                    marg(data(i,2),1) = data(i,1);
                end
                % Assign deterministic values from means of r.v.'s
                if isnan(data(i,3))
                    data(i,3) = marg(data(i,2),2);
                end
            end
            
            % Identifies r.v.'s, reliability parameters and deterministic parameters in probdata.marg
            Ix      = find(marg(:,1)>0);
            Ithetag = find(marg(:,1)==-1);
            Icg     = find(marg(:,1)==0);
            
            % Looks for reliability parameters in probdata.marg and transfers them to gfundata(lsf).thetag
            if ~isempty(Ithetag)
                for i=1:length(Ithetag)
                    if isempty(strmatch(margname(Ithetag(i)),dataname,'exact'))
                        thetagname = [ thetagname margname(Ithetag(i)) ];
                        thetag     = [ thetag; marg(Ithetag(i),2) ];
                    end
                end
            end
            
            % Looks for deterministic parameters in probdata.marg and transfers them to gfundata(lsf).cg
            if ~isempty(Icg)
                for i=1:length(Icg)
                    if isempty(strmatch(margname(Icg(i)),dataname,'exact'))
                        cgname = [ cgname margname(Icg(i)) ];
                        cg     = [ cg; marg(Icg(i),2) ];
                    end
                end
            end
            
            % Identifies r.v.'s, reliability parameters and deterministic parameters in femodel.data
            Ifex      = find(data(:,1)>0);
            Ifethetag = find(data(:,1)==-1);
            Ifecg     = find(data(:,1)==0);
            
            % Looks for reliability parameters in femodel.data and transfers them to gfundata(lsf).thetag
            % Sets Ithetag1, Ithetag2
            if ~isempty(Ifethetag)
                if isempty(thetag)
                    Ithetag1 = [];
                    if size(data,2)>3, thetagpara1 = []; end
                    if size(data,2)>4, thetagpara2 = []; end
                else
                    Ithetag1 = 1:length(thetag);
                    if size(data,2)>3, thetagpara1 = nan*ones(length(thetag),1); end
                    if size(data,2)>4, thetagpara2 = nan*ones(length(thetag),1); end
                end
                thetagname = [ thetagname dataname(Ifethetag) ];
                thetag     = [ thetag; data(Ifethetag,3) ];
                Ithetag2   = (length(Ithetag1)+1):(length(Ithetag1)+length(Ifethetag));
                if size(data,2)>3, thetagpara1 = [ thetagpara1; data(Ifethetag,4) ]; end
                if size(data,2)>4, thetagpara2 = [ thetagpara2; data(Ifethetag,5) ]; end
            else
                if isempty(thetag)
                    Ithetag1 = [];
                    if size(data,2)>3, thetagpara1 = []; end
                else
                    Ithetag1 = 1:length(thetag);
                    if size(data,2)>3, thetagpara1 = nan*ones(length(thetag),1); end
                end
                Ithetag2 = [];
                if size(data,2)>4, thetagpara2 = []; end
            end
            
            % Looks for deterministic parameters in femodel.data and transfer them into gfundata(lsf).cg
            % Sets Icg1, Icg2
            if ~isempty(Ifecg)
                if isempty(cg)
                    Icg1 = [];
                    if size(data,2)>3, cgpara1 = []; end
                    if size(data,2)>4, cgpara2 = []; end
                else
                    Icg1 = 1:length(cg);
                    if size(data,2)>3, cgpara1 = nan*ones(length(cg),1); end
                    if size(data,2)>4, cgpara2 = nan*ones(length(cg),1); end
                end
                cgname = [ cgname dataname(Ifecg) ];
                cg     = [ cg; data(Ifecg,3) ];
                Icg2   = (length(Icg1)+1):(length(Icg1)+length(Ifecg));
                if size(data,2)>3, cgpara1 = [ cgpara1; data(Ifecg,4) ]; end
                if size(data,2)>4, cgpara2 = [ cgpara2; data(Ifecg,5) ]; end
            else
                if isempty(cg)
                    Icg1 = [];
                    if size(data,2)>3, cgpara1 = []; end
                else
                    Icg1 = 1:length(cg);
                    if size(data,2)>3, cgpara1 = nan*ones(length(cg),1); end
                end
                Icg2 = [];
                if size(data,2)>4, cgpara2 = []; end
            end
            
            Ifexinmarg    = data(Ifex,2);
            Inonfexinmarg = (1:size(marg,1))';
            Inonfexinmarg([Ithetag; Icg; Ifexinmarg]) = [];
            
            Ix1 = 1:length(Inonfexinmarg);
            Ix2 = (length(Inonfexinmarg)+1):(length(Inonfexinmarg)+length(Ifexinmarg));
            
            if size(data,2)>3, xpara1 = [ nan*ones(length(Ix1),1); data(Ifex,4) ]; end
            if size(data,2)>4, xpara2 = [ nan*ones(length(Ix1),1); data(Ifex,5) ]; end
            
            margname    = margname([Inonfexinmarg; Ifexinmarg]);
            marg        = marg([Inonfexinmarg; Ifexinmarg],:);
            correlation = correlation([Inonfexinmarg; Ifexinmarg],[Inonfexinmarg; Ifexinmarg]);
            
        end % if ~isempty(femodel)
        
end % switch lower(gfundata(lsf).evaluator)


probdata.name        = margname;
probdata.marg        = marg;
probdata.correlation = correlation;
switch lower(gfundata(lsf).evaluator)
    case 'basic'
    otherwise
        if ~isempty(femodel)
            probdata.Ix1 = Ix1;
            probdata.Ix2 = Ix2;
            if size(data,2)>3, probdata.xpara1 = xpara1; end
            if size(data,2)>4, probdata.xpara2 = xpara2; end
        end % if ~isempty(femodel)
        
end


if isempty(thetag)
    if isfield(gfundata(lsf),'thetag')
        gfundata = rmfield(gfundata(lsf),{'thetag' 'thetagname'});
    end
else
    gfundata(lsf).thetagname = thetagname;
    gfundata(lsf).thetag     = thetag;
    switch lower(gfundata(lsf).evaluator)
        case 'basic'
        otherwise
            if ~isempty(femodel)
                gfundata(lsf).Ithetag1 = Ithetag1;
                gfundata(lsf).Ithetag2 = Ithetag2;
                if size(data,2)>3, gfundata(lsf).thetagpara1 = thetagpara1; end
                if size(data,2)>4, gfundata(lsf).thetagpara2 = thetagpara2; end
            end % if ~isempty(femodel)
    end
end


if isempty(cg)
    if isfield(gfundata(lsf),'cg')
        gfundata = rmfield(gfundata(lsf),{'cg' 'cgname'});
    end
else
    gfundata(lsf).cgname = cgname;
    gfundata(lsf).cg     = cg;
    switch lower(gfundata(lsf).evaluator)
        case 'basic'
        otherwise
            if ~isempty(femodel)
                gfundata(lsf).Icg1 = Icg1;
                gfundata(lsf).Icg2 = Icg2;
                if size(data,2)>3, gfundata(lsf).cgpara1 = cgpara1; end
                if size(data,2)>4, gfundata(lsf).cgpara2 = cgpara2; end
            end % if ~isempty(femodel)
    end
end

% Error check,
% thetaf sensitivity calculation is not implemented for some dsitribution function
% filters them out and throw an error message
if probdata.flag_sens == 1
    for i = 1:length(probdata.marg(:,1))
        type = probdata.marg(i,1);
        p4   = probdata.marg(i,8);
        % for these distribution types the sensitivities w.r.t. means, standard deviations, parameters and correlation coefficients
        % is not implemented
        if any(type == [20, 25, 30, 31])
            error(['For random variable number:',num2str(i), ' distribution type (20,30,31) the sensitivity calculation w.r.t. means, standard deviations, parameters and correlation coefficients is not yet implemented!'])
        end
        if all(type ~= [7, 51]) && not(p4 == 0 || isnan(p4))
            error(['For powered distribution, random variable number:',num2str(i), ', the sensitivity calculation w.r.t. means, standard deviations, parameters and correlation coefficients is not yet implemented!'])
        end
    end
end

% Error check,
% Nataf transformation is not implemented for some dsitribution function
% filters them out and throw an error message
for i = 1:length(probdata.marg(:,1))
    type  = probdata.marg(i,1);
    p4    = probdata.marg(i,8);
    corr_i = probdata.correlation(i,:);
    corr_i(i) = [];
    % check whether the random variable is correlated with some other
    if any(corr_i ~=0)
        % for these distribution types the Nataf transformation is not implemented
        if any(type == [20, 25, 30])
            error(['For random variable number: ',num2str(i), ' (has distribution type: ', num2str(type) ,') the Nataf transformation is not yet implemented!'])
        end
        % powered distributions
        if all(type ~= [7, 51]) && not(p4 == 0 || isnan(p4) || p4 == 1)
            error(['For powered distribution, random variable number: ',num2str(i), ', the Nataf transformation is not yet implemented!'])
        end
    end
end

analysisopt.already_updated = 1;