function  z = u_to_z(u,probdata)

% U_TO_Z    Transformation between u and z space
%
%   z = u_to_z(u,marg,varargin)
%
%   Function 'u_to_z' perform the transformation between standard
%   normal space and correlated normal space. 
%
%   Output: - z           = values of r.v.'s in correlated normal space
%   Input:  - u           = values of r.v.'s in standard normal space
%           - probdata

transf_type = probdata.transf_type;
marg        = probdata.marg;
nx          = size(u,2);
nrv         = size(marg,1);

x = zeros(nrv,nx);

switch transf_type
   
   case 3
      
      Lo = probdata.Lo;
      z = Lo * u;
      
   otherwise
      
end
