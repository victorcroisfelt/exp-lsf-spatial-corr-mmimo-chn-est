function [R] = functionExpLSF_ULA(M,theta,corrFactor,stdLSF)
%Generates the channel covariance matrix when considering spatially 
%correlated channels from the compound of the exponential spatial 
%correlation model and the consideration of large-scale fading (LSF) 
%variations over the array for a uniform linear array (ULA).
%
%This Matlab function is used in the paper:
%   
%Victor Croisfelt Rodrigues, Jose Carlos Marinello, and Taufik Abrao.
%"Exponential spatial correlation with large-scale fading variations in
%massive MIMO channel estimation". Trans Emerging Tel Tech. 2019;e3563.
%
%Download paper: https://doi.org/10.1002/ett.3563
%
%This is version 2.0 (Last edited: 04-09-2019)
%
%License: This code is licensed under the GPLv3 license. If you in any way
%use this code for research that results in publications, please reference 
%our original article as shown above.
%
%@Inputs:
% 	M: number of BS antennas.
%	theta: azimuthal angle.
%	corrFactor: correlation factor between antenna elements.
%   stdLSF: standard deviation of LSF variations over the array.
%
%@Outputs:
% 	R: M x M generated channel covariance matrix.
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%Generate realizations of LSF variations
LSFrealizations = randn(M,1);

%Put LSF realizations over the array as a diagonal matrix and in the power
%dimension
largeScaleFadingD = diag(10.^(stdLSF*LSFrealizations/20));

%Generate the covariance matrix with the exponential correlation model
%combined to LSF variations along the antenna elements
R = largeScaleFadingD*toeplitz((corrFactor*exp(1i*theta)).^(0:M-1))*largeScaleFadingD;