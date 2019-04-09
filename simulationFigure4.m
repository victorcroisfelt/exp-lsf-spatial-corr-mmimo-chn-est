%This Matlab script can be used to generate Figure 4 in the article:
%
%Victor Croisfelt Rodrigues, Jose Carlos Marinello, and Taufik Abrao.
%"Exponential Spatial Correlation with Large-Scale Fading Variations in
%Massive MIMO Channel Estimation". Trans Emerging Tel Tech. 2019;e3563.
%
%Download paper: https://doi.org/10.1002/ett.3563
%
%This is version 2.0 (Last edited: 04-09-2019)
%
%License: This code is licensed under the GPLv3 license. If you in any way
%use this code for research that results in publications, please reference 
%our original article as shown above.
%
%References:
%[1] Emil Bjornson, Jakob Hoydis and Luca Sanguinetti (2017), "Massive MIMO
%Networks: Spectral, Energy, and Hardware Efficiency", Foundations and
%Trends in Signal Processing: Vol. 11, No. 3-4, pp. 154-655. DOI: 10.1561/
%2000000093 (https://github.com/emilbjornson/massivemimobook).
%

%Initialization
close all;
clearvars;

%% Simulation parameters

%Choose the desired simulation subfigure:
%   simulation == 1: (a) r varying
%   simulation == 2: (b) sigma varying
%
simulation = 1;

%Number of BS antennas
M = 100;

if simulation == 1
    
    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = [0 1e-1:1e-1:1];
    
    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLSFrange = 0;
    
    %Extract the number of horizontal points
    nbrOfPoints = length(corrFactorRange);
    
elseif simulation == 2
    
    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = 0;
    
    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLSFrange = linspace(0,8,25);
    
    %Extract the number of horizontal points
    nbrOfPoints = length(stdLSFrange);
    
end

%Generate the range of nominal AoAs (Note: these variables define the Monte 
%Carlo average through the realizations of i.i.d. random  variables related
%to the randomness of the system; therefore, you must tuning these 
%parameters according to the desired accuracy)
thetaULAradians = linspace(-pi,+pi,10001);
thetaRadians = linspace(-pi,+pi,101);
varphiRadians = linspace(-pi/2,pi/2,101);

%Specify the effective SNR
SNRdB = 10;
SNR = 10.^(SNRdB/10);

%% Simulations

%Prepare to save the simulation results
NMSE_ULA = zeros(nbrOfPoints,length(thetaULAradians)-1);
NMSE_UPA = zeros(nbrOfPoints,length(thetaRadians)-1,length(varphiRadians)-1);

%Go through all horizontal points
for pts = 1:nbrOfPoints
    
    %Output simulation progress
    disp([num2str(pts) ' points out of ' num2str(nbrOfPoints)]);
    
    %Check the simulation choice and extract the current variables
    if simulation == 1
        
        corrFactor = corrFactorRange(pts);
        stdLSF = stdLSFrange;
        
    elseif simulation == 2
        
        corrFactor = corrFactorRange;
        stdLSF = stdLSFrange(pts);
        
    end
    
    %Go through all the azimuthal angles
    for az = 1:length(thetaULAradians)
        
        %Compute the spatial correlation matrix for the ULA
        R_ULA = functionExpLSF_ULA(M,thetaULAradians(az),corrFactor,stdLSF);
        
        %Compute the average squared norm of the estimation error
        C_ULA = R_ULA - SNR*R_ULA*((SNR*R_ULA+eye(M))\R_ULA);
        
        %Compute the NMSE
        NMSE_ULA(pts,az) = real(trace(C_ULA)/trace(R_ULA));
        
    end
    
    %Go through all azimuthal angles
    for az = 1:length(thetaRadians)-1
        
        %Go through all elevation angles
        for el = 1:length(varphiRadians)-1
            
            %Compute the spatial correlation matrix for the UPA
            R_UPA = functionExpLSF_UPA(M,thetaRadians(az),varphiRadians(el),corrFactor,stdLSF);
            
            %Compute the average squared norm of the estimation error
            C_UPA = R_UPA - SNR*R_UPA*((SNR*R_UPA+eye(M))\R_UPA);
            
            %Compute the NMSE
            NMSE_UPA(pts,az,el) = real(trace(C_UPA)/trace(R_UPA));
            
        end
        
    end
    
end

%Compute the NMSE for the uncorrelated fading case
NMSE_uncorrelated = 1/(SNR+1);

%% Plot simulation results

if simulation == 1
    
    figure;
    box on; hold on
    
    plot(corrFactorRange,NMSE_uncorrelated*ones(size(corrFactorRange)),'k:','LineWidth',1.5);
    plot(corrFactorRange,mean(NMSE_ULA,2),'--','LineWidth',1);
    plot(corrFactorRange,mean(mean(NMSE_UPA,2),3),'-','LineWidth',1);
    
    xlabel('Correlation factor ($r$)');
    ylabel('NMSE');
    
    legend('Bound: Uncorrelated Rayleigh Fading','ULA: Exponential w/ LSF variations','UPA: Exponential w/ LSF variations','Location','SouthWest');
       
    set(gca,'YScale','log');
     
elseif simulation == 2
    
    figure;
    box on; hold on
    
    plot(stdLSFrange,NMSE_uncorrelated*ones(size(stdLSFrange)),'k:','LineWidth',1.5);
    plot(stdLSFrange,mean(NMSE_ULA,2),'--','LineWidth',1);
    plot(stdLSFrange,mean(mean(NMSE_UPA,2),3),'-','LineWidth',1);
    
    xlabel('Standard deviation ($\sigma$)');
    ylabel('NMSE');
    
    legend('Bound: Uncorrelated Rayleigh Fading','ULA: Exponential w/ LSF variations','UPA: Exponential w/ LSF variations','Location','SouthWest');
    
    set(gca,'YScale','log');
    
    ylim([10^(-3) 10^(-1)])
      
end