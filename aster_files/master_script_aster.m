function master_script_aster(njobs0,njobs, ...
                             script_path, template_path, ...
                             mesh_generator, mesh_generator_arg, ...
                             work_parent, code_pathname, jobdir_name, job_name, ...
                             keep_optional_files,gawk_path_name)

% Perform a single call to the external code, here code_Aster FE code
% - External call to gmsh mesh generator
% - External call to code_Aster FE code
% - Postprocess code_Aster .resu file and generate gext_job file
                          
if njobs > 1
   disp('Error in master_script_aster.m: Switch analysisopt.multi_proc value to 0 in inputfile');
end
job = njobs0+1;

if strcmp(mesh_generator,'geo')
   
   % Create .geo file from template
   dos_cmd = [ gawk_path_name ' -f ' work_parent '/' jobdir_name num2str(job) '/replace_parameters_' num2str(job) '.awk < ' ...
                                     template_path '/template.geo > ' ...
                                     work_parent '/' jobdir_name num2str(job) '/' job_name '.geo' ];
   [cmd_flag,w] = dos(dos_cmd);
   
   % Execute gmsh                        
   dos_cmd = [ script_path '/call_gmsh.bat ' code_pathname ' ' work_parent '/' jobdir_name num2str(job) '/' job_name '.geo ' mesh_generator_arg ];
   [cmd_flag,w] = dos(dos_cmd);

elseif strcmp(mesh_generator,'msh')

elseif strcmp(mesh_generator,'dgibi')

elseif strcmp(mesh_generator,'mgib')

end

% Create .comm file from template
dos_cmd = [ gawk_path_name ' -f ' work_parent '/' jobdir_name num2str(job) '/replace_parameters_' num2str(job) '.awk < ' ...
                                  template_path '/template.comm > ' ...
                                  work_parent '/' jobdir_name num2str(job) '/' job_name '.comm' ];
[cmd_flag,w] = dos(dos_cmd);

% Create .export file from template
dos_cmd = [ gawk_path_name ' -f ' work_parent '/' jobdir_name num2str(job) '/replace_parameters_' num2str(job) '.awk < ' ...
                                  template_path '/template.export > ' ...
                                  work_parent '/' jobdir_name num2str(job) '/' job_name '.export' ];
[cmd_flag,w] = dos(dos_cmd);

% Execute code_Aster                          
dos_cmd = [ script_path '/call_aster.bat ' code_pathname ' ' work_parent '/' jobdir_name num2str(job) '/' job_name '.export' ];
[cmd_flag,w] = dos(dos_cmd);

% Extract required value from .resu file and store it in gext_job file                     
dos_cmd = [ gawk_path_name ' -f ' script_path '/extract_gext.awk < ' ...
                                  work_parent '/' jobdir_name num2str(job) '/' job_name '.resu > ' ...
                                  work_parent '/' jobdir_name num2str(job) '/gext_' num2str(job) ];
[cmd_flag,w] = dos(dos_cmd);

% Create gext file                     
dos_cmd = [ 'copy ' work_parent '/' jobdir_name num2str(job) '/gext_' num2str(job) ' ' work_parent '/gext' ];
dos_cmd(strfind(dos_cmd,'/'))='\';
[cmd_flag,w] = dos(dos_cmd);