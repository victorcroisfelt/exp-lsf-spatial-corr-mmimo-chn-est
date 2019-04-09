%This Matlab script can be used to generate Figure 8 in the article:
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

%Number of BS antennas
M = 100;

%Correlation factor (r) in the exponential correlation model
corrFactor = 0.5;

%Standard deviation [dB] of large-scale fading (LFS) variations over
%the array
stdLSF = 4;

%Define the number of statistical realizations (Note: this variable defines 
%the Monte Carlo average through the realizations of i.i.d. random 
%variables related to the randomness of the system; therefore, you must 
%tuning this parameter according to the desired accuracy)
nbrOfstats = 1e2;

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

%Define the effective SNR for the interfering UE (range)
SNR2dB = [10 0 -10];
SNR2 = 10.^(SNR2dB/10);

%% Simulation

%Prepare to store the simulation results
NMSE_ULA = zeros(length(thetaInterfererRadians),nbrOfstats,length(SNR2dB));
NMSE_UPA = zeros(length(thetaInterfererRadians),length(varphiInterfererRadians),nbrOfstats,length(SNR2dB));
NMSE_uncorr = zeros(2,length(SNR2dB));

%Go through all statistical realizations
for s = 1:nbrOfstats

    %Output simulation progress
    disp([num2str(s) ' stats out of ' num2str(nbrOfstats)]);

    %Compute the spatial correlation matrix of the desired UE for ULA 
    R1_ULA = functionExpLSF_ULA(M,desiredThetaUE1,corrFactor,stdLSF);

    %Compute the spatial correlation matrix of the desired UE for UPA 
    R1_UPA = functionExpLSF_UPA(M,desiredThetaUE1,desiredVarphiUE1,corrFactor,stdLSF);

    % Go through all azimuthal angles
    for az = 1:length(thetaInterfererRadians)

        %Compute the spatial correlation matrix of the interfering UE for ULA
        R2_ULA = functionExpLSF_ULA(M,thetaInterfererRadians(az),corrFactor,stdLSF); 

        %Go through all interfering SNRs
        for snr = 1:length(SNR2dB)

            %Compute the NMSE according (3.20) when having spatial correlation
            NMSE_ULA(az,s,snr) = 1 - SNR1*abs(trace(R1_ULA*((SNR1*R1_ULA+SNR2(snr)*R2_ULA+eye(M))\R1_ULA)))/trace(R1_ULA);

        end

        % Go through elevation angles
        for el = 1:length(varphiInterfererRadians)
        
            %Compute the spatial correlation matrix of the interfering UE for UPA
            R2_UPA = functionExpLSF_UPA(M,thetaInterfererRadians(az),varphiInterfererRadians(el),corrFactor,stdLSF); 
        
            %Go through all interfering SNRs
            for snr = 1:length(SNR2dB)
                
                %Compute the NMSE according (3.20) when having spatial correlation
                NMSE_UPA(az,el,s,snr) = 1 - SNR1*abs(trace(R1_UPA*((SNR1*R1_UPA+SNR2(snr)*R2_UPA+eye(M))\R1_UPA)))/trace(R1_UPA);
                
            end

        end
    
    end

end

%Go through all interfering SNRs
for snr = 1:length(SNR2dB)
    
    %Compute the NMSE according (3.20) when having uncorrelated fading
    NMSE_uncorr(:,snr) = 1 - SNR1/(SNR1+SNR2(snr)+1);
    
end

%% Plot the simulation results

%ULA

figure;
hold on; box on;

%Prepare to store y values
ycorr = zeros(length(thetaInterfererRadians),length(SNR2dB));

%Go through all interfering SNRs
for snr = 1:length(SNR2dB)

    %Extract the y values
    ycorr(:,snr) = mean(NMSE_ULA(:,:,snr),2);

end

%Plot the Uncorrelated Rayleigh Fading as reference
plot([thetaInterfererDegrees(1); thetaInterfererDegrees(end)],NMSE_uncorr(:,1),'k-','LineWidth',1.5)
plot([thetaInterfererDegrees(1); thetaInterfererDegrees(end)],NMSE_uncorr(:,2),'k--','LineWidth',1.5)
plot([thetaInterfererDegrees(1); thetaInterfererDegrees(end)],NMSE_uncorr(:,3),'k-.','LineWidth',1.5)

%Plot referece markers
plot(thetaInterfererDegrees(1),ycorr(1,1),'o-','LineWidth',1);
plot(thetaInterfererDegrees(1),ycorr(1,2),'o--','LineWidth',1);
plot(thetaInterfererDegrees(1),ycorr(1,3),'o-.','LineWidth',1);

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot the lines
plot(thetaInterfererDegrees,ycorr(:,1),'-','LineWidth',1);
plot(thetaInterfererDegrees,ycorr(:,2),'--','LineWidth',1);
plot(thetaInterfererDegrees,ycorr(:,3),'-.','LineWidth',1);

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot the missing markers
plot(thetaInterfererDegrees(20:20:end),ycorr((20:20:end),1),'o','LineWidth',1);
plot(thetaInterfererDegrees(20:20:end),ycorr((20:20:end),2),'o','LineWidth',1);
plot(thetaInterfererDegrees(20:20:end),ycorr((20:20:end),3),'o','LineWidth',1);

xlabel('Angle of interfering UE ($\theta$) [degree]');
ylabel('NMSE');

legend('Uncorrelated fading: same SNR','Uncorrelated fading: 10 dB weaker','Uncorrelated fading: 20 dB weaker','Correlated fading: same SNR','Correlated fading: 10 dB weaker','Correlated fading: 20 dB weaker','Location','SouthWest');

set(gca,'YScale','log');

xlim([-180 180]);
ylim([10^(-2) 10^(0)])

xticks(-180:45:180)


%UPA

figure;
hold on; box on; grid on;

%Compute the meshgrid space between the different angular dimensions
[vvarphiInterfererDegrees,tthetaInterfererDegrees] = meshgrid(varphiInterfererDegrees,thetaInterfererDegrees);

%Plot 3D-curve for UPA and same SNR
surf(tthetaInterfererDegrees,vvarphiInterfererDegrees,mean(NMSE_UPA(:,:,:,1),3),'FaceAlpha',0.8,'EdgeColor','b','LineWidth',2.0)

%Plot 3D-curve for UPA and 10 dB weaker
surf(tthetaInterfererDegrees,vvarphiInterfererDegrees,mean(NMSE_UPA(:,:,:,2),3),'FaceAlpha',0.8,'EdgeColor','b','LineWidth',2.0)

%Plot 3D-curve for UPA and 20 dB weaker
surf(tthetaInterfererDegrees,vvarphiInterfererDegrees,mean(NMSE_UPA(:,:,:,3),3),'FaceAlpha',0.8,'EdgeColor','b','LineWidth',2.0)

colormap(autumn);
shading interp

set(gca,'ZScale','log');

xlabel('$\theta$ [degree]');
ylabel('$\varphi$ [degree]');
zlabel('NMSE');

xlim([-180 180]);
ylim([-90 90]);
zlim([10^(-2) 10^(0)])

xticks(-180:45:180)
yticks(-90:15:90)

view(155,5)




