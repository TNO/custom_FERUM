function flag = purge_exec_dir(njobs0,njobs,femodel,client_OS)

% Purge the parent working directory

flag = 0;

% Computer platform parameters
script_path         = femodel.script_path;         if script_path(end) == '/', script_path = script_path(1:(end-1)); end
work_parent         = femodel.work_parent;         if work_parent(end) == '/', work_parent = work_parent(1:(end-1)); end
jobdir_name         = femodel.jobdir_name;

if strcmp(client_OS,'win')
   
   for job = (njobs0+1):(njobs0+njobs)
      
      cmd_flag = -1;
      n_trials = 0;
      while (cmd_flag ~= 0) & (n_trials < 20)
         dos_cmd = [work_parent '/' jobdir_name num2str(job)];
         dos_cmd(strfind(dos_cmd,'/'))='\';
         dos_cmd = ['rmdir /S /Q ' dos_cmd];
         [cmd_flag,w] = dos(dos_cmd);
         n_trials = n_trials+1;
         if cmd_flag ~= 0
            disp(['purge_exec_dir: error deleting content of directory ' work_parent '/' jobdir_name num2str(job) ' - retrying'])
            pause(1);
         end
      end
      
      if cmd_flag ~= 0
         disp(['purge_exec_dir: error deleting content of directory ' work_parent '/' jobdir_name num2str(job)])
         flag = -1;
         return
      end
      
   end
   
end