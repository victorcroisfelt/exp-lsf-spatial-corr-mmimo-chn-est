%This Matlab script can be used to generate Figure 5 in the article:
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

%Number of BS antennas (range)
Mrange = [4 16 100];

%Correlation factor (r) in the exponential correlation model
corrFactor = 0.5;

%Standard deviation [dB] of large-scale fading (LFS) variations over the 
%array
stdLSF = 4;

%Generate the range of nominal AoAs (Note: these variables define the Monte 
%Carlo average through the realizations of i.i.d. random  variables related
%to the randomness of the system; therefore, you must tuning these 
%parameters according to the desired accuracy)
thetaULAradians = linspace(-pi,+pi,10001);
thetaRadians = linspace(-pi,+pi,101);
varphiRadians = linspace(-pi/2,pi/2,101);

%Specify the effective SNR
SNRdB = (-10:20);
SNR = 10.^(SNRdB/10);

%% Simulations

%Prepare to save the simulation results
NMSE_ULA = zeros(length(SNR),length(thetaULAradians)-1,length(Mrange));
NMSE_UPA = zeros(length(SNR),length(thetaRadians)-1,length(varphiRadians)-1,length(Mrange));
NMSE_uncorr = zeros(length(SNR),1);

%Go through all SNR values (horizontal points)
for s = 1:length(SNR)
    
    %Output simulation progress
    disp([num2str(s) ' points out of ' num2str(length(SNR))]);
    
    %Compute the NMSE for the uncorrelated fading case
    NMSE_uncorr(s) = 1/(SNR(s)+1);
    
    %Go through all the azimuthal angles
    for az = 1:length(thetaULAradians)-1  
        
        %Compute the spatial correlation matrix for the ULA
        R_ULA = functionExpLSF_ULA(max(Mrange),thetaULAradians(az),corrFactor,stdLSF);
       
        %Go through all the BS antenna indexes
        for m = 1:length(Mrange)
            
            %Extract the correct R values
            Rm_ULA = R_ULA(1:Mrange(m),1:Mrange(m));
 
            %Compute the average squared norm of the estimation error
            C_ULA = Rm_ULA - SNR(s)*Rm_ULA*((SNR(s)*Rm_ULA+eye(Mrange(m)))\Rm_ULA);
            
            %Compute the NMSE
            NMSE_ULA(s,az,m) = real(trace(C_ULA)/trace(Rm_ULA));
            
        end
        
    end

    %Go through all azimuthal angles
    for az = 1:length(thetaRadians)-1

        %Go trough all elevation angles
        for el = 1:length(varphiRadians)-1

            %Compute the spatial correlation matrix for the UPA
            R_UPA = functionExpLSF_UPA(max(Mrange),thetaRadians(az),varphiRadians(el),corrFactor,stdLSF);

            %Go through all the BS antenna indexes
             for m = 1:length(Mrange)

                %Extract the R values
                Rm_UPA = R_UPA(1:Mrange(m),1:Mrange(m));

                %Compute the average squared norm of the estimation error
                C_UPA = Rm_UPA - SNR(s)*Rm_UPA*((SNR(s)*Rm_UPA+eye(Mrange(m)))\Rm_UPA);
            
                %Compute the NMSE
                NMSE_UPA(s,az,el,m) = real(trace(C_UPA)/trace(Rm_UPA));   

            end

        end

    end
    
end

%Compute the average values
avgNMSE_ULA_M1 = mean(NMSE_ULA(:,:,1),2);
avgNMSE_ULA_M16 = mean(NMSE_ULA(:,:,2),2);
avgNMSE_ULA_M100 = mean(NMSE_ULA(:,:,3),2);

avgNMSE_UPA_M1 = mean(mean(NMSE_UPA(:,:,:,1),2),3);
avgNMSE_UPA_M16 = mean(mean(NMSE_UPA(:,:,:,2),2),3);
avgNMSE_UPA_M100 = mean(mean(NMSE_UPA(:,:,:,3),2),3);

%% Plot simulation results

figure;
box on; hold on

%Plot referece markers for ULA
plot(SNRdB(1),avgNMSE_ULA_M1(1),'o--','LineWidth',1)
plot(SNRdB(1),avgNMSE_ULA_M16(1),'^--','LineWidth',1)
plot(SNRdB(1),avgNMSE_ULA_M100(1),'d--','LineWidth',1)

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot referece markers for UPA
plot(SNRdB(1),avgNMSE_UPA_M1(1),'o-','LineWidth',1)
plot(SNRdB(1),avgNMSE_UPA_M16(1),'^-','LineWidth',1)
plot(SNRdB(1),avgNMSE_UPA_M100(1),'d-','LineWidth',1)

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot the lines for ULA
plot(SNRdB,avgNMSE_ULA_M1,'--','LineWidth',1)
plot(SNRdB,avgNMSE_ULA_M16,'--','LineWidth',1)
plot(SNRdB,avgNMSE_ULA_M100,'--','LineWidth',1)

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot the lines for UPA
plot(SNRdB,avgNMSE_UPA_M1,'-','LineWidth',1)
plot(SNRdB,avgNMSE_UPA_M16,'-','LineWidth',1)
plot(SNRdB,avgNMSE_UPA_M100,'-','LineWidth',1)

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot the missing markers for ULA
plot(SNRdB(6:5:end),avgNMSE_ULA_M1(6:5:end),'o','LineWidth',1)
plot(SNRdB(6:5:end),avgNMSE_ULA_M16(6:5:end),'^','LineWidth',1)
plot(SNRdB(6:5:end),avgNMSE_ULA_M100(6:5:end),'d','LineWidth',1)

%Reset the colors
ax = gca;
ax.ColorOrderIndex = 1;

%Plot the missing markers for UPA
plot(SNRdB(6:5:end),avgNMSE_UPA_M1(6:5:end),'o','LineWidth',1)
plot(SNRdB(6:5:end),avgNMSE_UPA_M16(6:5:end),'^','LineWidth',1)
plot(SNRdB(6:5:end),avgNMSE_UPA_M100(6:5:end),'d','LineWidth',1)

%Plot the Uncorrelated Rayleigh Fading as reference
plot(SNRdB,NMSE_uncorr,'k:','LineWidth',1.5);

xlabel('Effective SNR');
ylabel('NMSE');

legend('ULA: $M = 4$','ULA: $M = 16$','ULA: $M = 100$','UPA: $M = 4$','UPA: $M = 16$','UPA: $M = 100$','Location','SouthWest');

set(gca,'YScale','log');

ylim([4e-3 10^(0)])
