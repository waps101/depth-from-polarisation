function [ height ] = HfPol( theta,diffuse,phi,l,mask,verbose,spec,azi,weight )
%HFPS Height from polarisation
% Inputs:
%   theta   - zenith angle from polarisation
%   diffuse - diffuse image (i.e. unpolarised appearance)
%   phi     - phase angle
%   (the above are all matrices size rows*cols)
%   l       - 3D column vector containing point light source
%   mask    - rows*cols binary foreground mask
%   verbose - if true, it will print out some details on the LS system
%   spec    - a binary spec mask (subset of mask)
%   azi     - priors on the azimuth, e.g. from a convexity prior
%   weight  - per-pixel weight for the azimuth prior
%
% Ouputs:
%   height  - rows*cols estimated height map
%
% William Smith 2016

% Todo: tidy up parameters so just pass an options parameter containing all
% optional params, with sensible defaults

% Note: this is a parameter - need to move this to options later
% Smoothing weight. lambda=0 => no smoothing.
lambda=0.1;

if nargin<6
    verbose=true;
end

rows = size(mask,1);
cols = size(mask,2);

f = cos(theta);

% Pad to avoid boundary problems
f = pad(f);
diffuse = pad(diffuse);
mask = pad(mask);
phi = pad(phi);
theta = pad(theta);
if nargin>=7
    spec = pad(spec);
else
    spec = zeros(size(mask));
end

if nargin==9
    weight = pad(weight);
    azi = pad(azi);
end

rows = rows+2;
cols = cols+2;

% Compute halfway vector
H = (l+[0 0 1]')./norm(l+[0 0 1]');
Hp = -H(1)./H(3);
Hq = -H(2)./H(3);

% Build lookup table relating x,y coordinate of valid pixels to index
% position in vectorised representation
count = 0;
indices = zeros(size(mask));
for row=1:rows
    for col=1:cols
        if mask(row,col)
            count=count+1;
            indices(row,col)=count;
        end
    end
end

% Create mask for 4-neighbours for Laplacian smoothing
h = [0 1 0; 1 0 1; 0 1 0];
mask4n = imfilter(1.*mask,h,'conv')==4;

% The number of usable pixels
npix = sum(mask(:));

if verbose
    disp(['Using ' num2str(npix) ' pixels']);
end

% Preallocate maximum required space
% This would be if all valid pixels had equations for all 8 neighbours for
% all possible pairs of images - it will be less in practice
i = zeros(npix*8*2+sum(mask4n(:))*4,1);
j = zeros(npix*8*2+sum(mask4n(:))*4,1);
s = zeros(npix*8*2+sum(mask4n(:))*4,1);

% Right hand side of linear system:
d = zeros(npix*2+sum(mask4n(:)),1);

NumEq = 0; % number of rows in matrix
k=0; % total number of non-zero entries in matrix

for row=1:rows
    for col=1:cols
        if mask(row,col)
            if mask4n(row,col) && lambda~=0
                % Add Laplacian smoothing term
                NumEq=NumEq+1;
                d(NumEq)=0;
                k=k+1;
                i(k)=NumEq; j(k)=indices(row,col); s(k)=-4*lambda;
                % Edge neighbours
                k=k+1;
                i(k)=NumEq; j(k)=indices(row,col-1); s(k)=1*lambda;
                k=k+1;
                i(k)=NumEq; j(k)=indices(row,col+1); s(k)=1*lambda;
                k=k+1;
                i(k)=NumEq; j(k)=indices(row-1,col); s(k)=1*lambda;
                k=k+1;
                i(k)=NumEq; j(k)=indices(row+1,col); s(k)=1*lambda;
            end
            if nargin==9
                if spec(row,col)
                    maxeq = 5;
                else
                    maxeq = 4;
                end
            else
                if spec(row,col)
                    maxeq = 3;
                else
                    maxeq = 2;
                end
            end
            for eq=1:maxeq
                
                % Start by computing p and q weights and RHSs of eqs
                if spec(row,col)
                    if eq==1 % Half vector constraint
                        RHS = Hp;
                        xval = 1;
                        yval = 0;
                    elseif eq==2 % Half vector constraint
                        RHS = Hq;
                        xval = 0;
                        yval = 1;
                    elseif eq==3 % Phase angle constraint
                        RHS = 0;
                        xval = cos(phi(row,col));
                        yval = -sin(phi(row,col));
                    elseif eq==4 % Azimuth prior from boundary (weighted)
                        RHS = -weight(row,col)*(sin(azi(row,col))*sin(theta(row,col)));
                        xval = weight(row,col)*cos(theta(row,col));
                        yval = 0;
                    elseif eq==5 % Azimuth prior from boundary (weighted)
                        RHS = -weight(row,col)*(cos(azi(row,col))*sin(theta(row,col)));
                        xval = 0;
                        yval = weight(row,col)*cos(theta(row,col));
                    end
                else
                    if eq==1 % Phase angle constraint
                        RHS = 0;
                        xval = cos(phi(row,col));
                        yval = -sin(phi(row,col));
                    elseif eq==2 % Ratio between DOP and Lambertian constraint
                        RHS = diffuse(row,col)/f(row,col) - l(3);
                        xval = -l(1);
                        yval = -l(2);
                    elseif eq==3 % Azimuth prior from boundary (weighted)
                        RHS = -weight(row,col)*(sin(azi(row,col))*sin(theta(row,col)));
                        xval = weight(row,col)*cos(theta(row,col));
                        yval = 0;
                    elseif eq==4 % Azimuth prior from boundary (weighted)
                        RHS = -weight(row,col)*(cos(azi(row,col))*sin(theta(row,col)));
                        xval = 0;
                        yval = weight(row,col)*cos(theta(row,col));
                    end
                end
                
                
                % Now decide which combination of neighbours are present
                % This determines which version of the numerical
                % approximation to the surface gradients will be used
                
                if mask(row,col-1) && mask(row,col+1)
                    % Both X neighbours present
                    if mask(row-1,col)
                        if mask(row+1,col)
                            if mask(row-1,col-1) && mask(row-1,col+1) && mask(row+1,col-1) && mask(row+1,col+1)
                                % All 8 neighbours present
                                NumEq=NumEq+1;
                                d(NumEq)=RHS;
                                % Edge neighbours
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-(4/12)*xval;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row,col+1); s(k)=(4/12)*xval;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row-1,col); s(k)=(4/12)*yval;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-(4/12)*yval;
                                % Corner neighbours
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row-1,col-1); s(k)=-(1/12)*xval+(1/12)*yval;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row-1,col+1); s(k)=(1/12)*xval+(1/12)*yval;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row+1,col-1); s(k)=-(1/12)*xval-(1/12)*yval;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row+1,col+1); s(k)=(1/12)*xval-(1/12)*yval;
                            else
                                % All 4 neighbours present
                                NumEq=NumEq+1;
                                d(NumEq)=RHS;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-xval/2;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row,col+1); s(k)=xval/2;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row-1,col); s(k)=yval/2;
                                k=k+1;
                                i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-yval/2;
                            end
                        else
                            % Both X, only forward in Y
                            NumEq=NumEq+1;
                            d(NumEq)=RHS;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-xval/2;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col+1); s(k)=xval/2;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row-1,col); s(k)=yval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col); s(k)=-yval;
                        end
                    elseif mask(row+1,col)
                        % Both X, only backward in Y
                        NumEq=NumEq+1;
                        d(NumEq)=RHS;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-xval/2;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col+1); s(k)=xval/2;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-yval;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col); s(k)=yval;
                    end
                elseif mask(row,col-1)
                    % Only backward in X
                    if mask(row-1,col)
                        if mask(row+1,col)
                            % Backward in X, both in Y
                            NumEq=NumEq+1;
                            d(NumEq)=RHS;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-xval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col); s(k)=xval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row-1,col); s(k)=yval/2;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-yval/2;
                        else
                            % Backward in X, only forward in Y
                            NumEq=NumEq+1;
                            d(NumEq)=RHS;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-xval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col); s(k)=xval-yval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row-1,col); s(k)=yval;
                        end
                    elseif mask(row+1,col)
                        % Backward in X, only backward in Y
                        NumEq=NumEq+1;
                        d(NumEq)=RHS;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-xval;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col); s(k)=xval+yval;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-yval;
                    end
                elseif mask(row,col+1)
                    % Only forward in X
                    if mask(row-1,col)
                        if mask(row+1,col)
                            % Forward in X, both in Y
                            NumEq=NumEq+1;
                            d(NumEq)=RHS;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col+1); s(k)=xval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col); s(k)=-xval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row-1,col); s(k)=yval/2;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-yval/2;
                        else
                            % Forward in X, only forward in Y
                            NumEq=NumEq+1;
                            d(NumEq)=RHS;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col+1); s(k)=xval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row,col); s(k)=-xval-yval;
                            k=k+1;
                            i(k)=NumEq; j(k)=indices(row-1,col); s(k)=yval;
                        end
                    elseif mask(row+1,col)
                        % Forward in X, only backward in Y
                        NumEq=NumEq+1;
                        d(NumEq)=RHS;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col+1); s(k)=xval;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row,col); s(k)=-xval+yval;
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-yval;
                    end
                end
            end
        end
        % Finished with a pixel
    end
end

if verbose
    disp(['System contains ' num2str(NumEq) ' linear equations with ' num2str(k) ' non-zero entries in C']);
end

i=i(1:k,1);
j=j(1:k,1);
s=s(1:k,1);
d=d(1:NumEq,1);

% Fix one pixel's height to zero to resolve unknown constant of integration
i(k+1)=NumEq+1;
j(k+1)=1;
s(k+1)=1;
d(NumEq+1)=0;

if sum(isnan(d))>0
    warning('d contains NaNs - output will be all NaN');
end

% Build matrix        +1 (for constraint on pixel 1)
C = sparse(i,j,s,NumEq+1,npix);
if verbose
    tic
end
z = C\d;
if verbose
    toc
end

% Copy vector of height values back into appropriate pixel positions
height=zeros(size(mask)).*NaN;
count=0;
for row=1:rows
    for col=1:cols
        if mask(row,col)
            count=count+1;
            height(row,col)=z(count);
        end
    end
end

% Unpad
height = height(2:rows-1,2:cols-1);

end

