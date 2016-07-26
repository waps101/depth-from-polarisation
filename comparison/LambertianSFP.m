function [ N,height ] = LambertianSFP( rho,phi,mask,n,s,albedo,Iun )
%LAMBERTIANSFP Shape-from-polarisation with Lambertian model
%   Inputs:
%      rho    - rows by cols matrix of DOP values
%      phi    - rows by cols matrix of phase angles
%      mask   - rows by cols binary foreground mask
%      n      - refractive index
%      s      - light source direction
%      albedo - diffuse albedo (can be scalar or per-pixel)
%      Iun    - unpolarised intensity
%
%   Outputs:
%      N      - rows by cols by 3 matrix containing surface normals
%      height - height map obtained by integrating N using lsqintegration
%
% This is a re-implementation of the technique used in:
%
% Mahmoud, A.H., El-Melegy, M.T., Farag, A.A.: Direct method for shape 
% recovery from polarization and shading. In: Proc. ICIP. (2012) 1769?1772
%
% William Smith
% 2016

% Invert degree of diffuse polarisation expression to compute zenith angle
temp = ((2.*rho + 2.*n.^2.*rho - 2.*n.^2 + n.^4 + rho.^2 + 4.*n.^2*rho.^2 - n.^4.*rho.^2 - 4.*n.^3.*rho.*(-(rho - 1).*(rho + 1)).^(1/2) + 1)./(n.^4.*rho.^2 + 2.*n.^4.*rho + n.^4 + 6.*n.^2.*rho.^2 + 4.*n.^2.*rho - 2.*n.^2 + rho.^2 + 2.*rho + 1)).^(1/2);
temp = min(real(temp),1);
theta = acos(temp);
Iun((Iun./albedo)>1)=albedo;
Iun = Iun./albedo;

% NOTE: In the paper, they find the intersection of the Lambertian
% constraint and the DOP constraint in terms of angles, I do it in terms of
% unit vectors. Once the two unit vectors have been found, the one that is
% closest to one of the ambiguous polarisation normals is chosen and the
% final result is the average of the chosen one and the closest
% polarisation normal

% Two possible ambiguous polarisation normals:
N1(:,:,1)=sin(phi).*sin(theta);
N1(:,:,2)=cos(phi).*sin(theta);
N1(:,:,3)=cos(theta);
N2(:,:,1)=sin(phi+pi).*sin(theta);
N2(:,:,2)=cos(phi+pi).*sin(theta);
N2(:,:,3)=cos(theta);

N = NaN(size(N1));
for row=1:size(mask,1)
    for col=1:size(mask,2)
        if mask(row,col)
            a = s(1);
            b = s(2);
            c = s(3)*cos(theta(row,col))-Iun(row,col);
            d = -sin(theta(row,col))^2;
            ny1 = -(b*c + a*(- d*a^2 - d*b^2 - c^2)^(1/2))/(a^2 + b^2);
            ny2 = -(b*c - a*(- d*a^2 - d*b^2 - c^2)^(1/2))/(a^2 + b^2);
            %if a~=0
            %    nx1 = -(c - (b*(b*c + a*(- d*a^2 - d*b^2 - c^2)^(1/2)))/(a^2 + b^2))/a;
            %    nx2 = -(c - (b*(b*c - a*(- d*a^2 - d*b^2 - c^2)^(1/2)))/(a^2 + b^2))/a;
            %else
                % Expression for nx1 and nx2 give Inf in this case so
                % compute from known norm of [nx ny]
                nx1 = sqrt(sin(theta(row,col))^2-ny1^2);
                nx2 = -sqrt(sin(theta(row,col))^2-ny1^2);
            %end
            % These are the two possible normals satisfying the Lambertian
            % + polarisation DOP constraint (intersection of two planes and
            % the unit norm constraint):
            n1 = [nx1 ny1 cos(theta(row,col))]';
            n2 = [nx2 ny2 cos(theta(row,col))]';
            n1best = max( dot(n1,squeeze(N1(row,col,:))), ...
                dot(n1,squeeze(N2(row,col,:))) );
            n2best = max( dot(n2,squeeze(N1(row,col,:))), ...
                dot(n2,squeeze(N2(row,col,:))) );
            if n1best>n2best
                na = n1;
            else
                na = n2;
            end
            if dot(na,squeeze(N1(row,col,:)))>dot(na,squeeze(N2(row,col,:)))
                nb = squeeze(N1(row,col,:));
            else
                nb = squeeze(N2(row,col,:));
            end
            n = na+nb;
            n(1:2) = n(1:2)./norm(n(1:2));
            n(1:2) = n(1:2).*sin(theta(row,col));
            n(3) = cos(theta(row,col));
            N(row,col,:)=real(n);
        end
    end
end

% Integrate normals into height map
height = lsqintegration( N,mask,false,[] );

end

