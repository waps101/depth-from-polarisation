function [ rho,phi,Iun ] = PolarisationImage( images,angles,mask,method )
%POLARISATIONIMAGE Decompose polarimetric images to polarisation image
%   Inputs:
%      images - rows by cols by nimages matrix of input images
%      angles - vector of polariser angles, length nimages
%      mask   - binary foreground mask (true for foreground)
%      method - either 'linear' or 'nonlinear'
%   Outputs:
%      rho    - rows by cols matrix containing degree of polarisation
%      phi    - rows by cols matrix containing phase angle
%      Iun    - rows by cols matrix containing unpolarised intensity
%
% William Smith 2016

if nargin<4
    % Default to linear method
    method='linear'
end

[rows,cols,nimages] = size(images);

if nargin>2
    % A mask has been provided so select only foreground pixels
    for i=1:nimages
        im = images(:,:,i);
        I(:,i)=im(mask);
    end
else
    % No mask - use all pixels
    I = reshape(images,rows*cols,nimages);
end

if strcmp(method,'nonlinear')
    % For each pixel solve a nonlinear optimisation problem 
    % Warning - slow (displays percentage progress)
    Iun = zeros(size(I,1),1);
    rho = zeros(size(I,1),1);
    phi = zeros(size(I,1),1);
    
    disp('      ');
    for i=1:size(I,1)
        [ Iun(i,1),rho(i,1),phi(i,1) ] = TRSfit( angles,I(i,:) );
        fprintf('\b\b\b\b\b\b%05.2f%%',i/size(I,1)*100);
    end
    phi = mod(phi,pi);
    
elseif strcmp(method,'linear')
    % Fast linear method - this seems to produce almost identical results
    % to the nonlinear method but is much faster
    A = [ones(nimages,1) cos(2.*angles') sin(2.*angles')];
    x = A\(I');
    x = x';
    Imax = x(:,1)+sqrt(x(:,2).^2+x(:,3).^2);
    Imin = x(:,1)-sqrt(x(:,2).^2+x(:,3).^2);
    Iun = (Imin+Imax)./2;
    rho = (Imax-Imin)./(Imax+Imin);
    phi = 0.5*atan2(x(:,3),x(:,2));
    phi = mod(phi,pi);
end

if nargin>2
    % We have a mask so need to reshape the estimated quantities to the
    % masked pixels only
    phi2 = zeros(rows,cols);
    phi2(mask) = phi;
    phi = phi2;
    rho2 = zeros(rows,cols);
    rho2(mask) = rho;
    rho = rho2;
    Iun2 = zeros(rows,cols);
    Iun2(mask) = Iun;
    Iun = Iun2;
else
    phi = reshape(phi,rows,cols);
    rho = reshape(rho,rows,cols);
    Iun = reshape(Iun,rows,cols);
end

end

