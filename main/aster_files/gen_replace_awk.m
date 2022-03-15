function [ source_path_names_str, destination_paths_str ] = gen_replace_awk(job,VA,expx,femodel,client_OS,source_path_names_str,destination_paths_str)

% Generate .awk replace files 

% Computer platform parameters
jobdir_name         = femodel.jobdir_name;
job_name            = femodel.job_name;

code_pathname       = femodel.code_pathname;       if code_pathname(end) == '/', code_pathname = code_pathname(1:(end-1)); end
script_path         = femodel.script_path;         if script_path(end) == '/', script_path = script_path(1:(end-1)); end
work_parent         = femodel.work_parent;         if work_parent(end) == '/', work_parent = work_parent(1:(end-1)); end
client_dir          = femodel.client_dir;          if client_dir(end) == '/', client_dir = client_dir(1:(end-1)); end

if strcmp(client_OS,'win')
   source_path_name      = [ work_parent '/' jobdir_name num2str(job) '/replace_parameters_' num2str(job) '.awk' ];   
   source_path_names_str = [];
   destination_paths_str = [];
end


% Generate .awk replace files 

if strcmp(client_OS,'win')

   fid_replace_awk = fopen(source_path_name,'wt');

   fprintf(fid_replace_awk, '{\n');

   fprintf(fid_replace_awk, 'gsub(/%%job%%/, "%d");\n', job);
   fprintf(fid_replace_awk, 'gsub(/%%jobdir_name%%/, "%s");\n', jobdir_name);
   fprintf(fid_replace_awk, 'gsub(/%%job_name%%/, "%s");\n', job_name);
   fprintf(fid_replace_awk, 'gsub(/%%code_pathname%%/, "%s");\n', code_pathname);
   fprintf(fid_replace_awk, 'gsub(/%%script_path%%/, "%s");\n', script_path);
   fprintf(fid_replace_awk, 'gsub(/%%work_parent%%/, "%s");\n', work_parent);

   for iexpx = 1:size(expx,1)
      fprintf(fid_replace_awk, 'gsub(/%%%s%%/, "%.10g");\n', VA(iexpx).name, expx(iexpx,1));
   end

   if isfield(femodel,'auxva') & femodel.auxva == 1
      auxva_string = auxva(VA,expx);
      nauxva       = length(auxva_string);
      for iauxva = 1:nauxva
         fprintf(fid_replace_awk, 'gsub(/%%auxva%d%%/, "%s");\n', iauxva, char(auxva_string(iauxva)) );
      end
   end

   fprintf(fid_replace_awk, 'print\n');
   fprintf(fid_replace_awk, '}\n');

   fclose(fid_replace_awk);

end
