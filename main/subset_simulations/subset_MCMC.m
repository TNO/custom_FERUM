function [ subtempu, I_eval, newG, ratio ] = subset_MCMC(simu,width,rand_generator)

% Generate a random walk from the current germ simu using the modified Metropolis-Hastings algorithm
% C.f. Au & Beck's paper


switch rand_generator
    case 0
        subtempu = rand(size(simu));
    otherwise
        % modern MATLAB has Mersenne Twister for generating pseudo
        % random numbers so we use that (FERUM's Twister throws an
        % error), `rand`'s default is the Twister algorithm
        subtempu = rand(size(simu));
end
subtempu = simu + ( (subtempu*width) - width/2 );


ratio = exp( -0.5*subtempu.^2 ) ./ exp( -0.5*simu.^2 );
I = find( ratio > 1 );
ratio(I) = 1;


switch rand_generator
    case 0
        test_unif = rand(size(subtempu)) >= ratio;
    otherwise
        % modern MATLAB has Mersenne Twister for generating pseudo
        % random numbers so we use that (FERUM's Twister throws an
        % error), `rand`'s default is the Twister algorithm
        test_unif = rand(size(subtempu)) >= ratio;
end


I = find( test_unif );
subtempu(I) = simu(I);


newG = ones(size(subtempu));
newG(I) = 0;

% Identify realizations of random vector to be evaluated
I_eval = find( any( ~test_unif , 1 ) );