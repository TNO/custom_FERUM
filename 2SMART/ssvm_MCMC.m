% Generate a random walk from the current germ Subgerm0 using the modified Metropolis-Hastings algorithm
% C.f. Au & Beck's paper

switch rand_generator
   case 0
      Subtemp = rand(size(Subgerm0));
   case 1
      Subtemp = twister(size(Subgerm0));
end
Subtemp = Subgerm0 + ( (Subtemp*width) - width/2 );

ratio  = exp( -0.5 * ((Subtemp-0)/1).^2 ) ./ exp( -0.5 * ((Subgerm0-0)/1).^2 );
I_McMc = find( ratio > 1 );
ratio(I_McMc) = 1;

switch rand_generator
   case 0
      test_unif = rand(size(Subtemp));
   case 1
      test_unif = twister(size(Subtemp));
end
  
I_McMc = find( test_unif >= ratio_amp_factor*ratio );
Subtemp(I_McMc) = Subgerm0(I_McMc);