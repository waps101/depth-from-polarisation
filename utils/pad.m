function [ X2 ] = pad( X )
%PAD Pad a matrix with an additional row/column at each edge

[rows,cols]=size(X);
X2 = zeros(rows+2,cols+2);
X2(2:rows+1,2:cols+1)=X(:,:);
if islogical(X)
    X2 = logical(X2);
end

end

