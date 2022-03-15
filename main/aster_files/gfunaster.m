function [ G, grad_g ] = gfunaster(lsf,x,grad_flag,probdata,analysisopt,gfundata,femodel,randomfield)

% Perform a single call to external FE code, here code_Aster
% - identify relevant environment variables
% - create all execution directories
% - generate all .awk replace files
% - external call to code_Aster FE code, through master_script_aster function
% - load gext file containing all gext_job values, assess G's corresponding limit-state values


global nfun

xname = probdata.name;
nx    = size(x,2);
nrv   = size(x,1);

% Random variables present in the FE model
Ix2   = probdata.Ix2;
nferv = length(Ix2);
% Other random variables
Ix1      = probdata.Ix1;
notherrv = length(Ix1);


for i = 1:nferv
   VA(i).name = char(probdata.name(Ix2(i)));
   if isfield(probdata,'xpara1')
      VA(i).xpara1 = probdata.xpara1(Ix2(i));
   end
   if isfield(probdata,'xpara2')
      VA(i).xpara2 = probdata.xpara2(Ix2(i));
   end
end


if isfield(gfundata(lsf),'thetag')
   
   thetagname    = gfundata(lsf).thetagname;
   thetag        = gfundata(lsf).thetag;
   Ithetag1      = gfundata(lsf).Ithetag1;
   Ithetag2      = gfundata(lsf).Ithetag2;
   nthetag       = size(gfundata(lsf).thetag,2);
   nthetagv      = size(gfundata(lsf).thetag,1);
   nfethetagv    = length(Ithetag2);
   notherthetagv = length(Ithetag1);
   for i = 1:nfethetagv
      VA(nferv+i).name = char(thetagname(Ithetag2(i)));
      if isfield(gfundata(lsf),'thetagpara1')
         VA(nferv+i).xpara1 = gfundata(lsf).thetagpara1(Ithetag2(i));
      end
      if isfield(gfundata(lsf),'thetagpara2')
         VA(nferv+i).xpara2 = gfundata(lsf).thetagpara2(Ithetag2(i));
      end
   end
   
else
   nthetag = 0; nfethetagv = 0; notherthetagv = 0;
   
end


if isfield(gfundata(lsf),'cg')
   
   cgname    = gfundata(lsf).cgname;
   cg        = gfundata(lsf).cg;
   Icg1      = gfundata(lsf).Icg1;
   Icg2      = gfundata(lsf).Icg2;
   ncg       = size(cg,2);
   ncgv      = size(cg,1);
   nfecgv    = length(Icg2);
   nothercgv = length(Icg1);
   for i = 1:nfecgv
      VA(nferv+nfethetagv+i).name = char(cgname(Icg2(i)));
      if isfield(gfundata(lsf),'cgpara1')
         VA(nferv+nfethetagv+i).xpara1 = gfundata(lsf).cgpara1(Icg2(i));
      end
      if isfield(gfundata(lsf),'cgpara2')
         VA(nferv+nfethetagv+i).xpara2 = gfundata(lsf).cgpara2(Icg2(i));
      end
   end
   
else
   
   ncg = 0; nfecgv = 0; nothercgv = 0;
   
end


if nx > 1

   EXPX = x(Ix2,:);
   if isfield(gfundata(lsf),'thetag')
      if ~isempty(Ithetag2)
         EXPX = [ EXPX; thetag(Ithetag2)*ones(1,nx) ];
      end
   end
   if isfield(gfundata(lsf),'cg')
      if ~isempty(Icg2)
	     EXPX = [ EXPX; cg(Icg2)*ones(1,nx) ];
      end
   end

elseif nthetag > 1

   % Not optimal: nthetag computations done, nfethetagv (<nthetag) required
   EXPX = x(Ix2)*ones(1,nthetag);
   if ~isempty(Ithetag2)
      EXPX = [ EXPX; thetag(Ithetag2,:) ];
   end
   if isfield(gfundata(lsf),'cg')
      if ~isempty(Icg2)
	     EXPX = [ EXPX; cg(Icg2)*ones(1,nthetag) ];
      end
   end

else

   EXPX = x(Ix2);
   if isfield(gfundata(lsf),'thetag')
      if ~isempty(Ithetag2)
	     EXPX = [ EXPX; thetag(Ithetag2) ];
	  end
   end
   if isfield(gfundata(lsf),'cg')
      if ~isempty(Icg2)
         EXPX = [ EXPX; cg(Icg2) ];
      end
   end

end


% Computer platform parameters

client_OS           = analysisopt.client_OS;

mesh_generator      = femodel.mesh_generator;
mesh_generator_arg  = femodel.mesh_generator_arg;  if isempty(mesh_generator_arg), mesh_generator_arg = 'nan'; end

jobdir_name         = femodel.jobdir_name;
job_name            = femodel.job_name;
keep_optional_files = femodel.keep_optional_files;

code_pathname       = femodel.code_pathname;       if code_pathname(end) == '/', code_pathname = code_pathname(1:(end-1)); end
script_path         = femodel.script_path;         if script_path(end) == '/', script_path = script_path(1:(end-1)); end
template_path       = femodel.template_path;       if template_path(end) == '/', template_path = template_path(1:(end-1)); end
work_parent         = femodel.work_parent;         if work_parent(end) == '/', work_parent = work_parent(1:(end-1)); end
client_dir          = femodel.client_dir;          if client_dir(end) == '/', client_dir = client_dir(1:(end-1)); end

gawk_path_name      = femodel.gawk_path_name;      if gawk_path_name(end) == '/', gawk_path_name = gawk_path_name(1:(end-1)); end


njobs0 = nfun;
if nx > 1
   njobs = nx;
elseif  nthetag > 1
   njobs = nthetag;
else
   njobs = 1;
end
block_size = analysisopt.block_size;


% Create all execution directories
flag = create_exec_dir(njobs0,njobs,femodel,client_OS);

% Generate all .awk replace files
source_path_names_str = '{';
destination_paths_str = '{';
for ijob = 1:njobs
   [ source_path_names_str, destination_paths_str ] = gen_replace_awk(njobs0+ijob,VA,EXPX(:,ijob),femodel,client_OS,source_path_names_str,destination_paths_str);
end
source_path_names_str = [ source_path_names_str ' }' ];
destination_paths_str = [ destination_paths_str ' }' ];

% External call to code_Aster
master_script_aster(njobs0,njobs, ...
                    script_path, template_path, ...
                    mesh_generator, mesh_generator_arg, ...
                    work_parent, code_pathname, jobdir_name, job_name, ...
                    keep_optional_files,gawk_path_name);

destination_gext_path = work_parent;

% Load gext file
eval([ 'gext = load(''' destination_gext_path '/gext'');' ]);
Gext = reshape(gext,1,njobs);

if grad_flag == 'yes'
   % Requires gext sensitivities (grad_gext) to be available from the FE code.
   % Needs not only grad_gext (output of the FE code used) but also differentiation
   % of g w.r.t. other r.v. appearing in g expression. Done through dgdq expressions
   eval([ 'grad_gext = load(''' destination_gext_path '/grad_gext'');' ]);
   Grad_gext = reshape(grad_gext,size(EXPX,1),njobs);
end

if keep_optional_files == 1
else
   % Purge the parent working directory
   flag = purge_exec_dir(njobs0,njobs,femodel,client_OS);
end


expression = gfundata(lsf).expression;


% Requires gext sensitivities (grad_gext) to be available from the FE code.
% Needs not only grad_gext (output of the FE code used) but also differentiation
% of g w.r.t. other r.v. appearing in g expression. Done through dgdq expressions
if grad_flag == 'yes'
   grad_expression = gfundata(lsf).dgdq;
end


if nx > 1

   G = zeros(1,nx);

   if isfield(gfundata(lsf),'thetag')
      thetagc = num2cell(thetag');
      eval(sprintf('[ %s ] = deal(thetagc{:});',[ char(thetagname)'; blanks(nthetagv)]));
   end

   if isfield(gfundata(lsf),'cg')
      cgc = num2cell(cg');
      eval(sprintf('[ %s ] = deal(cgc{:});',[ char(cgname)'; blanks(ncgv)]));
   end

   for i = 1:nx
       
      xic = num2cell(x(:,i)');
      eval(sprintf('[ %s ] = deal(xic{:});',[ char(xname)'; blanks(nrv)]));
      
      gext = Gext(i);
      G(i) = eval(expression);
      
   end

   grad_g = 0;

elseif nthetag > 1

   G = zeros(1,nthetag);
	
   xc = num2cell(x');
   eval(sprintf('[ %s ] = deal(xc{:});',[ char(xname)'; blanks(nrv)]));

   if isfield(gfundata(lsf),'cg')
      cgc = num2cell(cg');
      eval(sprintf('[ %s ] = deal(cgc{:});',[ char(cgname)'; blanks(ncgv)]));
   end

   for i = 1:nthetag
       
      thetagic = num2cell(thetag(:,i)');
      eval(sprintf('[ %s ] = deal(thetagic{:});',[ char(thetagname)'; blanks(nthetagv)]));
      
      gext = Gext(i);
      G(i) = eval(expression);
      
   end
	
   grad_g = 0;

else

   xc = num2cell(x');
   eval(sprintf('[ %s ] = deal(xc{:});',[ char(xname)'; blanks(nrv)]));

   if isfield(gfundata(lsf),'thetag')
      thetagc = num2cell(thetag');
      eval(sprintf('[ %s ] = deal(thetagc{:});',[ char(thetagname)'; blanks(nthetagv)]));
   end

   if isfield(gfundata(lsf),'cg')
      cgc = num2cell(cg');
      eval(sprintf('[ %s ] = deal(cgc{:});',[ char(cgname)'; blanks(ncgv)]));
   end

   gext = Gext;
   G = eval(expression);
	
   if grad_flag == 'yes'
      % Requires gext sensitivities (grad_gext) to be available from the FE code.
      % Needs not only grad_gext (output of the FE code used) but also differentiation
      % of g w.r.t. other r.v. appearing in g expression. Done through dgdq expressions
      for i = 1:nrv
         grad_gext = Grad_gext(i);
         grad_g(i,1) = eval(grad_expression(i));
      end
   else
      grad_g = 0;
   end

end

% Update number of calls to the limit state function
nfun = nfun + njobs;