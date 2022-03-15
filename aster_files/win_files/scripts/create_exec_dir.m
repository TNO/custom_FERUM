function flag = create_exec_dir(njobs0,njobs,femodel,client_OS)

% Create all the execution directories

flag = 0;

% Computer platform parameters
script_path         = femodel.script_path;         if script_path(end) == '/', script_path = script_path(1:(end-1)); end
work_parent         = femodel.work_parent;         if work_parent(end) == '/', work_parent = work_parent(1:(end-1)); end
jobdir_name         = femodel.jobdir_name;

if strcmp(client_OS,'win')
   
   for job = (njobs0+1):(njobs0+njobs)
      
      dos_cmd = ['mkdir ' work_parent '/' jobdir_name num2str(job)];
      dos_cmd(strfind(dos_cmd,'/'))='\';
      [cmd_flag,w] = dos(dos_cmd);
      if cmd_flag ~= 0
         disp(['create_exec_dir: error creating the directories ' work_parent '/' jobdir_name num2str(job)])
         flag = -1;
         return
      end

   end
   
end
