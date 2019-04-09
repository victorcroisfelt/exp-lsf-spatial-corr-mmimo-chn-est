%This Matlab script can be used to generate Figure 7 in the article:
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

%Number of BS antennas
M = 100;

%Correlation factor (r) in the exponential correlation model (range)
corrFactorRange = [0.5 1];

%Standard deviation [dB] of large-scale fading (LFS) variations over
%the array
stdLSF = 6;

%Define the number of statistical realizations (Note: this variable defines 
%the Monte Carlo average through the realizations of i.i.d. random 
%variables related to the randomness of the system; therefore, you must 
%tuning this parameter according to the desired accuracy)
nbrOfStats = 1e2;

%Specify the angles of the desired UE
desiredThetaUE1 = pi/6;
desiredVarphiUE1 = pi/6;

%Define the range of nominal angles of arrival (interfering UE)
thetaInterfererDegrees = -180:5:180;
thetaInterfererRadians = deg2rad(thetaInterfererDegrees);

varphiInterfererDegrees = -90:5:90;
varphiInterfererRadians = deg2rad(varphiInterfererDegrees);

%Define the effective SNR for the desired UE
SNR1dB = 10;
SNR1 = 10.^(SNR1dB/10);

%Define the effective SNR for the interfering UE
SNR2dB = 0;
SNR2 = 10.^(SNR2dB/10);

%% Simulation

%Prepare to store the simulation results
corrCoeff= zeros(length(thetaInterfererRadians),length(varphiInterfererRadians),nbrOfStats,length(corrFactorRange));

%Go through all different correlation scenarios
for scn = 1:length(corrFactorRange)

    %Output simulation progress
    disp([num2str(scn) ' scenarios out of ' num2str(length(corrFactorRange))]);

   %Go through all the statistical realizations
    for s = 1:nbrOfStats

    %Compute the spatial correlation matrix of the desired UE
    R1_UPA = functionExpLSF_UPA(M,desiredThetaUE1,desiredVarphiUE1,corrFactorRange(scn),stdLSF);
    
        %Go through all azimuthal angles
        for az = 1:length(thetaInterfererRadians)
            
            %Go through all elevation angles
            for el = 1:length(varphiInterfererRadians)
    
                %Compute the spatial correlation matrix for UPA of the interfering UE
                R2_UPA = functionExpLSF_UPA(M,thetaInterfererRadians(az),varphiInterfererRadians(el),corrFactorRange(scn),stdLSF);
                
                %Compute the denominator in (3.18)
                normalization = sqrt(SNR1*SNR2*abs(trace(R1_UPA*((SNR1*R1_UPA+SNR2*R2_UPA+eye(M))\R1_UPA)))*abs(trace(R2_UPA*((SNR1*R1_UPA+SNR2*R2_UPA+eye(M))\R2_UPA))));
                
                %Compute absolute value of antenna-averaged correlation coefficient in (3.18)
                corrCoeff(az,el,s,scn) = sqrt(SNR1*SNR2)*abs(trace(R1_UPA*((SNR1*R1_UPA+SNR2*R2_UPA+eye(M))\R2_UPA)))/normalization;
                
            end
            
        end
        
    end
    
end

%% Plot simulation results

figure;
hold on; box on; grid on;

%Compute the meshgrid space between the different angular dimensions
[vvarphiInterfererDegrees,tthetaInterfererDegrees] = meshgrid(varphiInterfererDegrees,thetaInterfererDegrees);

%Plot 3D-curve for UPA and r = 0.5 
surf(tthetaInterfererDegrees,vvarphiInterfererDegrees,(mean(corrCoeff(:,:,:,1),3)),'FaceAlpha', 0.8,'EdgeColor','b','LineWidth',2.0)

%Plot 3D-curve for UPA and r = 1
mesh(tthetaInterfererDegrees,vvarphiInterfererDegrees,(mean(corrCoeff(:,:,:,2),3)),'FaceAlpha', 0.8,'EdgeColor','b','LineWidth',2.0)

xlabel('$\theta$ [degree]');
ylabel('$\varphi$ [degree]');
zlabel('Antenna-averaged correlation coefficient');

xlim([-180 180]);
ylim([-90 90]);
zlim([0 0.2]);
