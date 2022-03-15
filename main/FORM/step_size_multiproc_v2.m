function step_size = step_size_multiproc_v2(lsf,G,grad_G,u,d,probdata,analysisopt,gfundata,femodel,randomfield)


global nfun

if isfield(analysisopt,'echo_flag')
   echo_flag = analysisopt.echo_flag;
else
   echo_flag = 1;
end

c = ( norm(u) / norm(grad_G) ) * 2 + 10;

merit = 0.5 * (norm(u))^2 + c * abs(G);

ntrial = 5;
Trial_step_size = 0.5.^[0:ntrial];

Trial_u = u * ones(1,ntrial+1) + d * Trial_step_size;

Trial_x = zeros(size(Trial_u));

for j = 1:(ntrial+1)
   trial_x = u_to_x(Trial_u(:,j),probdata);
	Trial_x(:,j) = trial_x;
end

switch analysisopt.multi_proc

	case 0

		trial_step_size = Trial_step_size(1);
        trial_u = Trial_u(:,1);
		trial_x = Trial_x(:,1);
		[ trial_G, ~ ] = gfun_v2(lsf,trial_x,'no',probdata,analysisopt,gfundata,femodel,randomfield);
		merit_new = 0.5 * (norm(trial_u))^2 + c * abs(trial_G);
		j = 1;
		while ( merit_new > merit ) && ( j < (ntrial+1) )
		   trial_step_size = Trial_step_size(1+j);
   		trial_u = Trial_u(:,1+j);
			trial_x = Trial_x(:,1+j);
			[ trial_G, ~ ] = gfun_v2(lsf,trial_x,'no',probdata,analysisopt,gfundata,femodel,randomfield);
			merit_new = 0.5 * (norm(trial_u))^2 + c * abs(trial_G);
			j = j + 1;
			if ( j == ntrial ) && ( merit_new > merit )
                if echo_flag
				   fprintf(1,'The step size has been reduced by a factor of 1/%d before continuing.',2^ntrial);
                end
			end
		end
		step_size = trial_step_size;

	case 1

		[ Trial_G, ~ ] = gfun(lsf,Trial_x,'no ',probdata,analysisopt,gfundata,femodel,randomfield);

		Merit_new = zeros(1,ntrial+1);
		for j = 1:(ntrial+1)
   		merit_new = 0.5 * (norm(Trial_u(:,j)))^2 + c * abs(Trial_G(j));
			Merit_new(j) = merit_new;
		end
		trial_step_size = Trial_step_size(1);
		merit_new = Merit_new(1);
		j = 1;
		while ( merit_new > merit ) && ( j < (ntrial+1) )
			trial_step_size = Trial_step_size(1+j);
			merit_new = Merit_new(1+j);
			j = j + 1;
			if ( j == (ntrial+1) ) && ( merit_new > merit )
                if echo_flag
                   fprintf(1,'The step size has been reduced by a factor of 1/%d.',2^ntrial);
                end
			end
		end
		step_size = trial_step_size;
end