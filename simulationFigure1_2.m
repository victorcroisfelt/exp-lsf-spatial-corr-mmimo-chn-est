%This Matlab script can be used to generate Figs. 1 and 2 in the article:
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
    corrFactorRange = [0 0.25 0.5 0.75 0.99 1];
    
    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLFSrange = 0;
    
    %Extract the number of different evaluated scenarios
    nbrOfScenarios = length(corrFactorRange);
    
elseif simulation == 2
    
    %Correlation factor (r) in the exponential correlation model (range)
    corrFactorRange = 0;
    
    %Standard deviation [dB] of large-scale fading (LFS) variations over
    %the array (range)
    stdLFSrange = [0 2 4 6 8];
    
    %Extract the number of different evaluated scenarios
    nbrOfScenarios = length(stdLFSrange);
    
end

%Generate the range of nominal AoAs (Note: these variables define the Monte 
%Carlo average through the realizations of i.i.d. random  variables related
%to the randomness of the system; therefore, you must tuning these 
%parameters according to the desired accuracy)
thetaULAradians = linspace(-pi,+pi,10001); %aiming the average equality
thetaRadians = linspace(-pi,pi,101); %azimuth
varphiRadians = linspace(-pi/2,pi/2,101); %elevation

%% Simulation

%Prepare to save simulation results
eigenvalues_ULA = zeros(M,length(thetaULAradians)-1,nbrOfScenarios);
eigenvalues_UPA = zeros(M,length(thetaRadians)-1,length(varphiRadians)-1,nbrOfScenarios);

%Go through all different scenarios
for scn = 1:nbrOfScenarios
    
    %Display simulation progress
    disp([num2str(scn) ' scenarios out of ' num2str(nbrOfScenarios)]);
    
    if simulation == 1
        
        %Extract the current correlation factor value
        corrFactor = corrFactorRange(scn);
        stdLSF = stdLFSrange; %just repeat the stdLSF
        
    elseif simulation == 2
        
        %Extract the current LSF standard deviation value
        stdLSF = stdLFSrange(scn);
        corrFactor = corrFactorRange; %just repeat the correlation factor
        
    end
    
    %Go through all azimuthal angles
    for az = 1:length(thetaULAradians)-1
        
        %Generate covariance matrix with the exponential correlation model with
        %fading variations along the array for the uniform linear array (ULA)
        R_ULA = functionExpLSF_ULA(M,thetaULAradians(az),corrFactor,stdLSF);
        
        %Save simulation results for ULA
        eigenvalues_ULA(:,az,scn) = sort(eig(R_ULA),'descend');
        
    end
    
    %Go through all azimuthal angles
    for az = 1:length(thetaRadians)-1
        
        %Go through all elevation angles
        for el = 1:length(varphiRadians)-1
            
            %Generate covariance matrix with the exponential correlation model with
            %fading variations along the array for the uniform planar array (UPA)
            R_UPA = functionExpLSF_UPA(M,thetaRadians(az),varphiRadians(el),corrFactor,stdLSF);
            
            %Save simulation results for UPA
            eigenvalues_UPA(:,az,el,scn) = sort(eig(R_UPA),'descend');
            
        end
        
    end
    
end

%% Plot simulation results

if simulation == 1
    
    %ULA
    
    figure;
    hold on; box on;
    
    %Prepare to store y values
    y = zeros(M,nbrOfScenarios);
    
    %Go through all different scenarios
    for scn = 1:nbrOfScenarios
        
        %Extract the y values
        y(:,scn) = real(mean(eigenvalues_ULA(:,:,scn),2));
        
    end
    
    %Plot referece markers
    plot(1,y(1,1),'h--','LineWidth',1);
    plot(1,y(1,2),'*--','LineWidth',1);
    plot(1,y(1,3),'s--','LineWidth',1);
    plot(1,y(1,4),'d--','LineWidth',1);
    plot(1,y(1,5),'^--','LineWidth',1);
    plot(1,y(1,6),'p--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the lines
    plot(1:M,y(:,1),'--','LineWidth',1);
    plot(1:M,y(:,2),'--','LineWidth',1);
    plot(1:M,y(:,3),'--','LineWidth',1);
    plot(1:M,y(:,4),'--','LineWidth',1);
    plot(1:M,y(:,5),'--','LineWidth',1);
    plot(1:M,y(:,6),'--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the missing markers
    plot(10:10:M,y(10:10:M,1),'h','LineWidth',1);
    plot(10:10:M,y(10:10:M,2),'*','LineWidth',1);
    plot(10:10:M,y(10:10:M,3),'s','LineWidth',1);
    plot(10:10:M,y(10:10:M,4),'d','LineWidth',1);
    plot(10:10:M,y(10:10:M,5),'^','LineWidth',1);
    plot(10:10:M,y(10:10:M,6),'p','LineWidth',1);
    
    %Plot the uncorrelated reference line
    plot(1:M,ones(M,1),'k:','LineWidth',1.5);
    text(55,1.25,'Uncorrelated Rayleigh Fading')
   
    xlabel('Eigenvalue number (decreasing order)');
    ylabel('Average normalized eigenvalue');
    
    legend('$r = 0$','$r = 0.25$','$r = 0.5$','$r = 0.75$','$r = 0.99$','r = 1','Location','NorthEast');
    
    set(gca,'Yscale','log');

    ylim([1e-3 1e3]);
    
    clear y
    
    %UPA
    
    figure;
    hold on; box on;
    
    %Prepare to store y values
    y = zeros(M,nbrOfScenarios);
    
    %Go through all different scenarios
    for scn = 1:nbrOfScenarios
        
        %Extract the y values
        y(:,scn) = real(mean(mean(eigenvalues_UPA(:,:,:,scn),2),3));
        
    end
    
    %Plot referece markers
    plot(1,y(1,1),'h--','LineWidth',1);
    plot(1,y(1,2),'*--','LineWidth',1);
    plot(1,y(1,3),'s--','LineWidth',1);
    plot(1,y(1,4),'d--','LineWidth',1);
    plot(1,y(1,5),'^--','LineWidth',1);
    plot(1,y(1,6),'p--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the lines
    plot(1:M,y(:,1),'--','LineWidth',1);
    plot(1:M,y(:,2),'--','LineWidth',1);
    plot(1:M,y(:,3),'--','LineWidth',1);
    plot(1:M,y(:,4),'--','LineWidth',1);
    plot(1:M,y(:,5),'--','LineWidth',1);
    plot(1:M,y(:,6),'--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the missing markers
    plot(10:10:M,y(10:10:M,1),'h','LineWidth',1);
    plot(10:10:M,y(10:10:M,2),'*','LineWidth',1);
    plot(10:10:M,y(10:10:M,3),'s','LineWidth',1);
    plot(10:10:M,y(10:10:M,4),'d','LineWidth',1);
    plot(10:10:M,y(10:10:M,5),'^','LineWidth',1);
    plot(10:10:M,y(10:10:M,6),'p','LineWidth',1);
    
    %Plot the uncorrelated reference line
    plot(1:M,ones(M,1),'k:','LineWidth',1.5);
    text(55,1.25,'Uncorrelated Rayleigh Fading')
   
    xlabel('Eigenvalue number (decreasing order)');
    ylabel('Average normalized eigenvalue');
    
    legend('$r = 0$','$r = 0.25$','$r = 0.5$','$r = 0.75$','$r = 0.99$','r = 1','Location','NorthEast');
    
    set(gca,'Yscale','log');

    ylim([1e-3 1e3]);
    
elseif simulation == 2
    
    %ULA
    
    figure;
    hold on; box on;
    
    %Prepare to store y values
    y = zeros(M,nbrOfScenarios);
    
    %Go through all different scenarios
    for scn = 1:nbrOfScenarios
        
        %Extract the y values
        y(:,scn) = real(mean(eigenvalues_ULA(:,:,scn),2));
        
    end
    
    %Plot referece markers
    plot(1,y(1,1),'h--','LineWidth',1);
    plot(1,y(1,2),'*--','LineWidth',1);
    plot(1,y(1,3),'s--','LineWidth',1);
    plot(1,y(1,4),'d--','LineWidth',1);
    plot(1,y(1,5),'^--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the lines
    plot(1:M,y(:,1),'--','LineWidth',1);
    plot(1:M,y(:,2),'--','LineWidth',1);
    plot(1:M,y(:,3),'--','LineWidth',1);
    plot(1:M,y(:,4),'--','LineWidth',1);
    plot(1:M,y(:,5),'--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the missing markers
    plot(10:10:M,y(10:10:M,1),'h','LineWidth',1);
    plot(10:10:M,y(10:10:M,2),'*','LineWidth',1);
    plot(10:10:M,y(10:10:M,3),'s','LineWidth',1);
    plot(10:10:M,y(10:10:M,4),'d','LineWidth',1);
    plot(10:10:M,y(10:10:M,5),'^','LineWidth',1);
    
    %Plot the uncorrelated reference line
    plot(1:M,ones(M,1),'k:','LineWidth',1.5);
    text(55,1.25,'Uncorrelated Rayleigh Fading')
   
    xlabel('Eigenvalue number (decreasing order)');
    ylabel('Average normalized eigenvalue');
    
    legend('$\sigma = 0$ dB','$\sigma = 2$ dB','$\sigma = 4$ dB','$\sigma = 6$ dB','$\sigma = 8$ dB','Location','NorthEast');    
    
    set(gca,'Yscale','log');

    ylim([1e-3 1e3]);
    
    %UPA
    
    figure;
    hold on; box on;
    
    %Prepare to store y values
    y = zeros(M,nbrOfScenarios);
    
    %Go through all different scenarios
    for scn = 1:nbrOfScenarios
        
        %Extract the y values
        y(:,scn) = real(mean(mean(eigenvalues_UPA(:,:,:,scn),2),3));
        
    end
    
    %Plot referece markers
    plot(1,y(1,1),'h--','LineWidth',1);
    plot(1,y(1,2),'*--','LineWidth',1);
    plot(1,y(1,3),'s--','LineWidth',1);
    plot(1,y(1,4),'d--','LineWidth',1);
    plot(1,y(1,5),'^--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the lines
    plot(1:M,y(:,1),'--','LineWidth',1);
    plot(1:M,y(:,2),'--','LineWidth',1);
    plot(1:M,y(:,3),'--','LineWidth',1);
    plot(1:M,y(:,4),'--','LineWidth',1);
    plot(1:M,y(:,5),'--','LineWidth',1);
    
    %Reset the colors
    ax = gca;
    ax.ColorOrderIndex = 1;
    
    %Plot the missing markers
    plot(10:10:M,y(10:10:M,1),'h','LineWidth',1);
    plot(10:10:M,y(10:10:M,2),'*','LineWidth',1);
    plot(10:10:M,y(10:10:M,3),'s','LineWidth',1);
    plot(10:10:M,y(10:10:M,4),'d','LineWidth',1);
    plot(10:10:M,y(10:10:M,5),'^','LineWidth',1);
    
    %Plot the uncorrelated reference line
    plot(1:M,ones(M,1),'k:','LineWidth',1.5);
    text(55,1.25,'Uncorrelated Rayleigh Fading')

    xlabel('Eigenvalue number (decreasing order)');
    ylabel('Average normalized eigenvalue');
    
    legend('$\sigma = 0$ dB','$\sigma = 2$ dB','$\sigma = 4$ dB','$\sigma = 6$ dB','$\sigma = 8$ dB','Location','NorthEast');
    
    set(gca,'Yscale','log');

    ylim([1e-3 1e3]);
    
end