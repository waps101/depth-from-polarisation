function [ l,T,B ] = findLight( theta,phi,diffuse,mask,ldim,l )
%FINDLIGHT Estimate illumination from polarisation data
% Inputs:
%   theta   - zenith angle estimates from degree of polarisation
%   phi     - phase angle from polarisation image
%   diffuse - unpolarised intensity
%   mask    - binary foreground mask
%   ldim    - Set to 3, 4 or 9 depending on whether solving for point
%             source, order 1 or order 2 spherical harmonic illumination
%   l       - if ldim=3 and the direction of l is known (but not the
%             intensity) then pass in the known direction and only the
%             intensity will be optimised
% Outputs:
%   l       - 3D, 4D or 9D lighting coefficient vector
%   T       - Ambiguous transformation, T*l is also a solution
%   B       - The basis that is consistent with l, B*T is also a possible
%             basis
%
% William Smith 2016

if nargin<5
    ldim = 3;
end

%% Set up some required values

i = diffuse(mask);
theta = theta(mask);
phi = phi(mask);
N(:,1)=sin(phi).*sin(theta);
N(:,2)=cos(phi).*sin(theta);
N(:,3)=cos(theta);

%% Set up appropriate (ambiguous) basis vectors and transformation matrices depending on dimensionality of lighting

if ldim==3
    B = N;
    T = [-1 0 0; 0 -1 0; 0 0 1];
elseif ldim==4
    B = N;
    B(:,4) = 1;
    T = diag([-1 -1 1 1]);
elseif ldim==9
    B = [ones(size(N,1),1) N 3.*N(:,3).^2-1 N(:,1).*N(:,2) N(:,1).*N(:,3) N(:,2).*N(:,3) N(:,1).^2-N(:,2).^2];
    T = diag([1 -1 -1 1 1 1 -1 -1 1]);
end

%% Initialise optimisation

converged = false;
tau = 1e-9;
niter = 0;
maxiter = 100;

if nargin==6
    % In this case, we already know the light source direction and just
    % want to optimise for the intensity
    while ~converged
        idx = (B*l-i).^2>(B*T*l-i).^2;
        P = B;
        P(idx,:)=P(idx,:)*T;
        lnew=l.*((P*l)\i);
        if norm(l-lnew)<tau
                converged=true;
        end
        l=lnew;
        niter = niter+1;
        if niter>maxiter
            break
        end
    end
else
    % Random initialisation
    l = randn(ldim,1);
    if ldim==4
        % This is treating the order 1 model as ambient plus point source
        % and therefore assuming ambient should be positive
        l(4)=abs(l(4));
    end
    
    %% Perform alternating assignment/optimisation
    
    while ~converged
        % indicator function
        idx = (B*l-i).^2>(B*T*l-i).^2;
        P = B;
        P(idx,:)=P(idx,:)*T;
        if ldim==3
            % Unconstrained linear least squares
            lnew = P\i;
        elseif ldim==4
            % Order 0 term is constrained to be positive
            lnew = lsqlin(P,i,[],[],[],[],[-inf -inf -inf 0],[inf inf inf inf]);
        elseif ldim==9
            lnew = P\i;
        end
        if norm(l-lnew)<tau
            converged=true;
        end
        l = lnew;
        niter = niter+1;
        if niter>maxiter
            break
        end
    end
end
B = P;

disp(['Number of iterations = ' num2str(niter)]);
disp(['Residual = ' num2str(norm(P*l-i).^2)]);

end

