function [R] = functionExpLSF_UPA(M,theta,varphi,corrFactor,stdLSF)
%Generates the channel covariance matrix when considering spatially 
%correlated channels from the compound of the exponential spatial 
%correlation model and the consideration of large-scale fading (LSF) 
%variations over the array for a uniform planar array (UPA).
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
%   varphi: elevation angle.
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

%Check for the special case of M = 1
if M == 1
  
    %Determine the spatial correlation matrix of the horizontal ULA
    R = functionExpLSF_ULA(M,theta,corrFactor,stdLSF);   
    
else

    %Compute the number of antennas at the horizontal and vertical 
    %dimensions, where M must be a perfect square number
    M = sqrt(M);

    %Determine the spatial correlation matrix of the horizontal ULA
    Rh = functionExpLSF_ULA(M,theta,corrFactor,stdLSF);

    %Determine the spatial correlation matrix of the vertical ULA
    Rv = functionExpLSF_ULA(M,varphi,corrFactor,stdLSF);

    %Compute the generated R covariance matrix by using the Kronecker model
    R = kron(Rh,Rv);
    
end
