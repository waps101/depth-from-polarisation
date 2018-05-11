function [ theta_s ] = rho_spec( rho,n )
%RHO_S Compute zenith angle from degree of polarisation for specular pix
%   There is no closed form inversion for rho_s so this functions uses a
%   look up table and interpolates

theta = 0:0.01:pi/2;

rho_s = (2.*sin(theta).^2.*cos(theta).*sqrt(n.^2-sin(theta).^2))./(n.^2-sin(theta).^2-n.^2.*sin(theta).^2+2.*sin(theta).^4);
%figure; plot(theta,rho_s)
maxpos = find(rho_s==max(rho_s));

theta = theta(1:maxpos);
rho_s = rho_s(1:maxpos);

theta_s = interp1(rho_s,theta,rho);

end

