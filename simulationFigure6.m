%This Matlab script can be used to generate Figure 6 in the article:
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

%Number of BS antennas (range)
Mrange = [1 4 16 100];

if simulation == 1

    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = [0.5 1];

    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLSFrange = 0;

    %Extract the number of different evaluated scenarios
    nbrOfScenarios = length(corrFactorRange);
    
    %Specify the number of statistical realizations
    nbrOfStats = 1;

elseif simulation == 2

    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = 0.5;

    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLSFrange = [2 6];

    %Extract the number of different evaluated scenarios
    nbrOfScenarios = length(stdLSFrange);
      
    %Define the number of statistical realizations (Note: this variable 
    %defines the Monte Carlo average through the realizations of i.i.d. 
    %random variables related to the randomness of the system; therefore, 
    %you must tuning this parameter according to the desired accuracy)
    nbrOfStats = 1e2;

end

%Specify the angles of the desired UE
desiredThetaUE1 = pi/6;

%Define the range of nominal angles of arrival
thetaInterfererDegrees = -180:1:180;
thetaInterfererRadians = deg2rad(thetaInterfererDegrees);

%Define the effective SNR for the desired UE
SNR1dB = 10;
SNR1 = 10.^(SNR1dB/10);

%Define the effective SNR for the interfering UE
SNR2dB = 0;
SNR2 = 10.^(SNR2dB/10);

%% Simulation

%Prepare to store the simulation results
corrCoeff = zeros(length(thetaInterfererRadians),nbrOfStats,nbrOfScenarios,length(Mrange));

%Go through all different correlation scenarios
for scn = 1:nbrOfScenarios

    %Output simulation progress
    disp([num2str(scn) ' scenarios out of ' num2str(nbrOfScenarios)]);

    %Check the simulation choice and extract the current variables
    if simulation == 1

        corrFactor = corrFactorRange(scn);
        stdLSF = stdLSFrange;

    elseif simulation == 2

        corrFactor = corrFactorRange;
        stdLSF = stdLSFrange(scn);

    end
    
    %Go through all the statistical realizations
    for s = 1:nbrOfStats

    %Compute the spatial correlation matrix of the desired UE for ULA
    R1_ULA = functionExpLSF_ULA(max(Mrange),desiredThetaUE1,corrFactor,stdLSF);

        %Go through all azimuthal angles
        for az = 1:length(thetaInterfererRadians)

            %Compute the spatial correlation matrix for ULA of the infering 
            %UE
            R2_ULA = functionExpLSF_ULA(max(Mrange),thetaInterfererRadians(az),corrFactor,stdLSF);

            %Go through all BS antenna indexes
            for m = 1:length(Mrange)

                %Extract correlation matrices of the specified dimension
                R1m = R1_ULA(1:Mrange(m),1:Mrange(m));
                R2m = R2_ULA(1:Mrange(m),1:Mrange(m));
                
                %Compute the denominator in (3.18)
                normalization = sqrt(SNR1*SNR2*abs(trace(R1m*((SNR1*R1m+SNR2*R2m+eye(Mrange(m)))\R1m)))*abs(trace(R2m*((SNR1*R1m+SNR2*R2m+eye(Mrange(m)))\R2m))));
                
                %Compute absolute value of antenna-averaged correlation 
                %coefficient in (3.18)
                corrCoeff(az,s,scn,m) = sqrt(SNR1*SNR2)*abs(trace(R1m*((SNR1*R1m+SNR2*R2m+eye(Mrange(m)))\R2m)))/normalization;

            end
            
        end
        
    end
    
end
    
%% Plot simulation results
  
if simulation == 1

    figure;
    hold on; box on;

    %Plot curves for ULA and r = 0.5 
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,1),2),'-','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,2),2),'--','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,3),2),'-.','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,4),2),':','LineWidth',1)
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot curves for ULA and r = 1
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,1),2),'-','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,2),2),'--','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,3),2),'-.','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,4),2),':','LineWidth',1)

    %Plot the Uncorrelated Rayleigh Fading as reference
    plot(thetaInterfererDegrees,ones(size(thetaInterfererDegrees)),'k:','LineWidth',1.5);
    text(-170,1.024,'Uncorrelated Rayleigh Fading')

    xlabel('Angle of the interfering UE [degree]');
    ylabel('Antenna-averaged correlation coefficient');

    legend('$M=1$','$M=4$','$M=16$','$M=100$','Location','SouthWest');

    xlim([-180 180]);
    ylim([0 1.1]);

elseif simulation == 2

    figure;
    hold on; box on;

    %Plot curves for ULA and sigma = 2 dB 
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,1),2),'-','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,2),2),'--','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,3),2),'-.','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,1,4),2),':','LineWidth',1)
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot curves for ULA and sigma = 6 dB 
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,1),2),'-','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,2),2),'--','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,3),2),'-.','LineWidth',1)
    plot(thetaInterfererDegrees,mean(corrCoeff(:,:,2,4),2),':','LineWidth',1)    
   
    plot(thetaInterfererDegrees,ones(size(thetaInterfererDegrees)),'k:','LineWidth',1.5);
    text(-170,1.024,'Uncorrelated Rayleigh Fading')

    %Plot the Uncorrelated Rayleigh Fading as reference
    xlabel('Angle of the interfering UE [degree]');
    ylabel('Antenna-averaged correlation coefficient');

    legend('$M=1$','$M=4$','$M=16$','$M=100$','Location','Best');

    xlim([-180 180]);
    ylim([0 1.1]);

end

