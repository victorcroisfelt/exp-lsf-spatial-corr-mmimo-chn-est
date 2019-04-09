%This Matlab script can be used to generate Figure 3 in the article:
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

%Define the range of the number of BS antennas
Mrange = (1:27).^2;
Mmax = max(Mrange); % extract the max value

if simulation == 1

    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = [0.5 0.75];

    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLSFrange = 0;
    
    %Extract the number of different evaluated scenarios
    nbrOfScenarios = length(corrFactorRange);
    
    %Define number of statistical realizations
    nbrOfStats = 1;
    
elseif simulation == 2
    
    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = 0;

    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLSFrange = [4 6];

    %Extract the number of different evaluated scenarios
    nbrOfScenarios = length(stdLSFrange);
    
    %Define the number of statistical realizations (Note: this variable 
    %defines the Monte Carlo average through the realizations of i.i.d. 
    %random variables related to the randomness of the system; therefore, 
    %you must tuning this parameter according to the desired accuracy)
    nbrOfStats = 1e2; 
    
end

%Define the position of desired UE
desiredTheta = pi/6;
desiredVarphi = pi/6;

%% Simulation

%Prepare to store simulation results  
variance_ULA = zeros(length(Mrange),nbrOfStats,nbrOfScenarios); 
variance_UPA = zeros(length(Mrange),nbrOfStats,nbrOfScenarios);
    
%Go through all different scenarios
for scn = 1:nbrOfScenarios
   
    %Output simulation progress
    disp([num2str(scn) ' scenarios out of ' num2str(nbrOfScenarios)]); 
    
    %Extract the current values
    if simulation == 1
       
        corrFactor = corrFactorRange(scn);
        stdLSF = stdLSFrange;
        
    elseif simulation == 2
        
        corrFactor = corrFactorRange;
        stdLSF = stdLSFrange(scn);
        
    end
    
    %Prepare to store the several covariance matrix realizations
    R_ULA = zeros(Mmax,Mmax,nbrOfStats);
    R_UPA = zeros(Mmax,Mmax,nbrOfStats);
    
    %Go through all statistics realizations
    for s = 1:nbrOfStats
    
        %Generates the spatial correlation matrices
        R_ULA(:,:,s) = functionExpLSF_ULA(Mmax,desiredTheta,corrFactor,stdLSF);
        R_UPA(:,:,s) = functionExpLSF_UPA(Mmax,desiredTheta,desiredVarphi,corrFactor,stdLSF);
        
    end
    
    %Go through all BS antenna values
    for m = 1:length(Mrange)
        
        %Extract values 
        Rm_ULA = R_ULA(1:Mrange(m),1:Mrange(m),:);
        Rm_UPA = R_UPA(1:Mrange(m),1:Mrange(m),:);
        
        %Go through all statistics realizations
        for s = 1:nbrOfStats
            
            %Compute the variance of the channel hardening
            variance_ULA(m,s,scn) = real(trace(Rm_ULA(:,:,s)*Rm_ULA(:,:,s))/(trace(Rm_ULA(:,:,s))).^2);
            variance_UPA(m,s,scn) = real(trace(Rm_UPA(:,:,s)*Rm_UPA(:,:,s))/(trace(Rm_UPA(:,:,s))).^2);
                
        end
       
    end
   
end

%% Plot simulation results

if simulation == 1
    
    figure;
    hold on; box on;
    
    %Plot the variance of Uncorrelated Rayleigh Fading as reference
    plot(Mrange,1./Mrange,'k:','LineWidth',1.5);
    
    %Plot results for ULA
    plot(Mrange,mean(variance_ULA(:,:,1),2),'o--','LineWidth',1);
    plot(Mrange,mean(variance_ULA(:,:,2),2),'^--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot results for UPA
    plot(Mrange,mean(variance_UPA(:,:,1),2),'o-','LineWidth',1);
    plot(Mrange,mean(variance_UPA(:,:,2),2),'^-','LineWidth',1);
    
    xlabel('Number of BS antennas (M)');
    ylabel('$v_{jjk}$');

    legend('Bound: Uncorrelated Rayleigh Fading','ULA: Correlated, $r = 0.5$','ULA: Correlated, $r = 0.75$','UPA: Correlated, $r = 0.5$','UPA: Correlated, $r = 0.75$','Location','NorthEast');
    
    xlim([0 200])
    ylim([0 0.25])
    
elseif simulation == 2
   
    figure;
    hold on; box on;
    
    %Plot the variance of Uncorrelated Rayleigh Fading as reference
    plot(Mrange,1./Mrange,'k:','LineWidth',1.5);
    
    %Plot the results for ULA
    plot(Mrange,mean(variance_ULA(:,:,1),2),'o--','LineWidth',1);
    plot(Mrange,mean(variance_ULA(:,:,2),2),'^--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the results for UPA
    plot(Mrange,mean(variance_UPA(:,:,1),2),'o-','LineWidth',1);
    plot(Mrange,mean(variance_UPA(:,:,2),2),'^-','LineWidth',1);
    
    xlabel('Number of BS antennas ($M$)');
    ylabel('Channel hardening variance metric ($v_{k}$)');

    legend('Bound: Uncorrelated Rayleigh Fading','ULA: Correlated, $\sigma = 4$ dB','ULA: Correlated, $\sigma = 6$ dB','UPA: Correlated, $\sigma = 4$ dB','UPA: Correlated, $\sigma = 6$ dB','Location','NorthEast');
    
    xlim([0 200])
    ylim([0 0.25])
    
end