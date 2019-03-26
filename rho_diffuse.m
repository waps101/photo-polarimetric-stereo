function [ theta ] = rho_diffuse( rho,n )
%RHO_DIFFUSE Compute zenith angle from degree of polarisation for diffuse pix
%   Uses closed form inversion of DOP expression

temp = ((2.*rho + 2.*n.^2.*rho - 2.*n.^2 + n.^4 + rho.^2 + 4.*n.^2*rho.^2 - n.^4.*rho.^2 - 4.*n.^3.*rho.*(-(rho - 1).*(rho + 1)).^(1/2) + 1)./(n.^4.*rho.^2 + 2.*n.^4.*rho + n.^4 + 6.*n.^2.*rho.^2 + 4.*n.^2.*rho - 2.*n.^2 + rho.^2 + 2.*rho + 1)).^(1/2);
% To avoid complex result in case of numerical issues:
temp = min(real(temp),1);
theta = acos(temp);

end

