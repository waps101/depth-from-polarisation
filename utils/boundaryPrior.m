function [ azi,Bdist ] = boundaryPrior( mask,weight )
%BOUNDARYPRIOR Use foreground mask to compute convexity prior
%   Input:
%      mask   - rows by cols binary foreground mask
%      weight - determines how fast the per-pixel weight falls off from the
%               boundary (default: 5)
%   Output:
%      azi    - rows by cols azimuth angle prior
%      Bdist  - weight associated with each azimuth estimate

if nargin<2
    weight=5;
end

% Find boundary of mask
B = bwperim(mask);
[bcol,brow]=find(B);
[row,col]=meshgrid(1:size(mask,2),1:size(mask,1));
% Compute distance to closest point on boundary
[~,D]=knnsearch([brow bcol],[row(mask) col(mask)]);
Bdist = zeros(size(mask)).*NaN;
% Transform to range 0..1
Bdist(mask)=-D;
Bdist(mask) = Bdist(mask)-min(Bdist(mask));
Bdist = Bdist./max(Bdist(mask));
Bdist = Bdist.^weight;
mask2 = mask;
azi = zeros(size(mask));
% Repeatedly erode boundary to propagate into interior
while sum(mask2(:))>0
    % Find boundary of current mask
    B = bwboundaries(mask2,8);
    B{1}(end+1,:)=B{1}(1,:);
    for i=1:size(B{1},1)-1
        azi(B{1}(i,1),B{1}(i,2))=atan2(B{1}(i+1,2)-B{1}(i,2),B{1}(i+1,1)-B{1}(i,1));
    end
    % Remove boundary pixels from mask, i.e. erode mask
    mask2(bwperim(mask2,8))=0;
end
% Transform azimuth angles to vectors
dx = cos(azi);
dy = sin(azi);
dx(~mask)=0;
dy(~mask)=0;
% Smooth vector field
g=fspecial('Gaussian',[9 9],2);
dx = imfilter(dx,g);
dy = imfilter(dy,g);
% Transform back to azimuth angle
azi = atan2(dy,dx);
azi = mod(-azi+pi-pi/2,2*pi);

end

