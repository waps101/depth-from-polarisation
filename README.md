# Linear depth estimation from an uncalibrated, monocular polarisation image

This is a Matlab implementation of our ECCV 2016 paper "Linear depth estimation from an uncalibrated, monocular polarisation image". It also includes an implementation of polarimetric image decomposition (linear and nonlinear optimisation), two comparison shape-from-polarisation methods, a simple least squares surface integration method (with foreground mask) and a basic method for pixel-wise specular labelling.

Note: I am in the process of cleaning up the code and adding to the repository. I will update the list of what has been uploaded as I go along. Content included so far:

1. Comparison methods
2. Least squares integrator
3. Polarimetric image decomposition
4. ... check back soon for more!

I will add documentation and demo scripts as I upload the code.

## Polarimetric image decomposition

The first thing you need to do is convert your captured image into a 3-channel polarisation image. The function that does this is PolarisationImage.m. Inputs are:

1. images - 3D array containing captured images of size rows by cols by nimages
2. angles - vector of length nimages containing polariser angles (I use a coordinate system where the polariser angle is measured from the upward vertical axis, increasing in a clockwise direction if viewed looking into the camera lens)
3. (optional) mask - binary foreground mask of size rows by cols
4. (optional) method - either 'linear' or 'nonlinear', default: linear

It returns rho (degree of polarisation), phi (phase angle) and Iun (unpolarised intensity).

Sample call:

```matlab
[ rho,phi,Iun ] = PolarisationImages( images,angles,mask,'linear' );
figure; imagesc(rho); colorbar
figure; imagesc(phi); colorbar
figure; imshow(Iun)
```
