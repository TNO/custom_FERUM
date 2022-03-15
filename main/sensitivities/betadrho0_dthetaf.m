function f = betadrho0_dthetaf(drho0_dthetaf,dF_dthetaf,PHI2,F,dPHI2_dZ,dZ_dthetaf,dPHI2_drho0,detJ,WIP)

%f = sum( sum( ( dF_dthetaf .* PHI2 + F .* ( dPHI2_dZ .* dZ_dthetaf + dPHI2_drho0 * drho0_dthetaf ) ) * detJ .* WIP ) );
f = sum( sum( ( dF_dthetaf .* PHI2 + F .* ( dPHI2_drho0 * drho0_dthetaf ) ) * detJ .* WIP ) );
