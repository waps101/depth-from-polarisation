function [ Iun,rho,phi ] = TRSfit( angles,I )
%TRSFIT Nonlinear least squares optimisation to fit sinusoid
%   Inputs:
%      angles - vector of polarising filter angles
%      I      - vector of measured intensities
%   Outputs:
%      Iun, rho, phi - scalar values containing polarisation image params
%
% William Smith 2016

% Initialise
b0(1)=mean(I);
b0(2)=sqrt(mean((I-b0(1)).^2)).*sqrt(2);
b0(3)=0;

options = optimoptions('lsqnonlin', 'Algorithm','trust-region-reflective','Display','off');
b=lsqnonlin(@(b) TRSfitobj(b,angles,I),b0,[0 0 -pi],[inf inf pi],options);

if b(3)<0
    b(3)=b(3)+pi;
end
Iun = b(1);
Imax = Iun+b(2);
Imin = Iun-b(2);
rho = (Imax-Imin)/(Imax+Imin);
phi = b(3);

end

function errs = TRSfitobj(b,angles,I)

I2 = b(1)+b(2).*cos(2.*angles-2*b(3));
errs = I-I2;

end