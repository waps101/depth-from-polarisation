% Load raw images, mask and specular mask
load sampleData.mat

% Estimate polarisation image from captured images
[ rho_est,phi_est,Iun_est ] = PolarisationImage( images,angles,mask,'linear' );

% Assume refractive index = 1.5
n = 1.5;

% Estimate light source direction from diffuse pixels (note that you might
% get a convex/concave flip)
%[ s,T,B ] = findLight( theta_est,phi_est,Iun_est,mask&~spec,3 );
% Or use known direction and estimate albedo
s = [2 0 7]';
[ s,T,B ] = findLight( theta_est,phi_est,Iun_est,mask&~spec,3,s );

% Compute angles, taking into account different model for specular pixels
theta_est_combined = rho_diffuse(rho_est,n);
theta_s = rho_spec(rho_est(spec),n);
theta_est_combined(spec)=theta_s;
phi_est_combined = phi_est;
phi_est_combined(spec)=mod(phi_est(spec)+pi/2,pi);

% Compute boundary prior azimuth angles and weight
[ azi,Bdist ] = boundaryPrior( mask );

% Run linear height from polarisation
[ height ] = HfPol( theta_est_combined,min(1,Iun_est),phi_est_combined,s,mask,false,spec );

% Visualise
figure;
surf(height,'EdgeColor','none','FaceColor',[0 0 1],'FaceLighting','gouraud','AmbientStrength',0,'DiffuseStrength',1); axis equal; light