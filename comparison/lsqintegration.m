function [ height ] = lsqintegration( N,mask,verbose,guidez )
%LSQINTEGRATION Least squares surface integration with optional guide z
%   Input:
%      N      - rows by cols by 3 matrix containing surface normals
%      mask   - rows by cols binary foreground mask
%      guidez - rows by cols guide depth map
%   Output:
%      z      - estimated depth map
%
% William Smith 2016

if nargin<3
    verbose=true;
end

[rows,cols,~]=size(N);

% Pad to avoid boundary problems
N2 = zeros(rows+2,cols+2,3);
N2(2:rows+1,2:cols+1,:)=N;
N = N2;
clear N2

P = -N(:,:,1)./N(:,:,3);
Q = -N(:,:,2)./N(:,:,3);

mask2 = zeros(rows+2,cols+2);
mask2(2:rows+1,2:cols+1)=mask;
mask = mask2;
clear mask2

rows = rows+2;
cols = cols+2;

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

% Create mask for smoothed central difference filter
h = [1 0 1; 1 0 1; 1 0 1];
maskSCDx = imfilter(1.*mask,h,'conv')==6;
maskSCDy = imfilter(1.*mask,h','conv')==6;
SCDx = (1/12).*[-1 0 1; -4 0 4; -1 0 1];
SCDy = flipud(SCDx');

% Create mask for SavGol gradient filter
h = ones(5,5);
h(:,3)=0;
maskSGx = imfilter(1.*mask,h,'conv')==20;
maskSGy = imfilter(1.*mask,h','conv')==20;
[~, ~, ~, ~, SG] = SavGol(3,5);
SGx = SG(:,:,2)';
SGy = flipud(SGx');

% Create mask for 4-neighbours for Laplacian smoothing
h = [0 1 0; 1 0 1; 0 1 0];
mask4n = imfilter(1.*mask,h,'conv')==4;
% Smoothing weight. lambda=0 => no smoothing.
lambda=0.2;

% The number of usable pixels
npix = sum(mask(:));

if verbose
    disp(['Using ' num2str(npix) ' pixels']);
end

% Preallocate maximum required space
% This would be if all valid pixels had equations for all 8 neighbours for
% all possible pairs of images - it will be less in practice
i = zeros(npix*20*2,1);
j = zeros(npix*20*2,1);
s = zeros(npix*20*2,1);

% Right hand side of linear system:
d = zeros(npix*2,1);

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
            % X equations
%             if maskSGx(row,col) && 1==0
%                 NumEq=NumEq+1;
%                 d(NumEq)=P(row,col);
%                 for a=1:5
%                     for b=1:5
%                         if SGx(a,b)~=0
%                         k=k+1;
%                         i(k)=NumEq; j(k)=indices(row+a-3,col+b-3); s(k)=SGx(a,b);
%                         end
%                     end
%                 end
%             elseif
            if maskSCDx(row,col)
                NumEq=NumEq+1;
                d(NumEq)=P(row,col);
                for a=1:3
                    for b=1:3
                        if SCDx(a,b)~=0
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+a-2,col+b-2); s(k)=SCDx(a,b);
                        end
                    end
                end
            elseif mask(row,col-1)
                if mask(row,col+1)
                    % Central difference
                    NumEq=NumEq+1;
                    d(NumEq)=P(row,col);
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col+1); s(k)=0.5;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-0.5;
                else
                    % Backward difference
                    NumEq=NumEq+1;
                    d(NumEq)=P(row,col);
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col-1); s(k)=-1;
                end
            elseif mask(row,col+1)
                % Forward difference
                NumEq=NumEq+1;
                d(NumEq)=P(row,col);
                k=k+1;
                i(k)=NumEq; j(k)=indices(row,col); s(k)=-1;
                k=k+1;
                i(k)=NumEq; j(k)=indices(row,col+1); s(k)=1;
            end
            
            % Y equations
%             if maskSGy(row,col) && 1==0
%                 NumEq=NumEq+1;
%                 d(NumEq)=Q(row,col);
%                 for a=1:5
%                     for b=1:5
%                         if SGy(a,b)~=0
%                         k=k+1;
%                         i(k)=NumEq; j(k)=indices(row+a-3,col+b-3); s(k)=SGy(a,b);
%                         end
%                     end
%                 end
%             elseif 
            if maskSCDy(row,col)
                NumEq=NumEq+1;
                d(NumEq)=Q(row,col);
                for a=1:3
                    for b=1:3
                        if SCDy(a,b)~=0
                        k=k+1;
                        i(k)=NumEq; j(k)=indices(row+a-2,col+b-2); s(k)=SCDy(a,b);
                        end
                    end
                end
            elseif mask(row-1,col)
                if mask(row+1,col)
                    % Central difference
                    NumEq=NumEq+1;
                    d(NumEq)=Q(row,col);
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row-1,col); s(k)=0.5;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-0.5;
                else
                    % Backward difference
                    NumEq=NumEq+1;
                    d(NumEq)=Q(row,col);
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row,col); s(k)=-1;
                    k=k+1;
                    i(k)=NumEq; j(k)=indices(row-1,col); s(k)=1;
                end
            elseif mask(row+1,col)
                % Forward difference
                NumEq=NumEq+1;
                d(NumEq)=Q(row,col);
                k=k+1;
                i(k)=NumEq; j(k)=indices(row,col); s(k)=1;
                k=k+1;
                i(k)=NumEq; j(k)=indices(row+1,col); s(k)=-1;
            end
        end
    end
end

if verbose
    disp(['System contains ' num2str(NumEq) ' linear equations with ' num2str(k) ' non-zero entries in C']);
end

i=i(1:k,1);
j=j(1:k,1);
s=s(1:k,1);
d=d(1:NumEq,1);

i(k+1)=NumEq+1;
j(k+1)=1;
s(k+1)=1;
d(NumEq+1)=0;

% Build matrix        +1 (for constraint on pixel 1)
C = sparse(i,j,s,NumEq+1,npix);
if verbose
    tic
end
%sum(isnan(C(:)))
%sum(isnan(d(:)))
%[c,A] = qr(C,d,0);
%sum(isnan(A(:)))
%sum(isnan(c(:)))

%z = A\c;
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
        
