
close all;
clear;
L = 16;
K = 10;
Mrange = 10:10:100;
Mmax = max(Mrange);
Reuse_f = 4;
nbrOfSetups = 10;
nbrOfRealizations = 10;

%% Propagation parameters
B = 20e6;
p_UL = 100;
p_DL = 100;
noiseFigure = 7;
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;
tau_c = 200;
accuracy = 2;
SD = 10;
sumSE_MR = zeros(length(Mrange),length(Reuse_f),nbrOfSetups);
sumSE_RZF = zeros(length(Mrange),length(Reuse_f),nbrOfSetups);
sumSE_MMMSE = zeros(length(Mrange),length(Reuse_f),nbrOfSetups);

%% Go through all setups
for n = 1:nbrOfSetups
    f=4;
    disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
    [R,channelGaindB] = functionExampleSetup(L,K,Mmax,accuracy,SD);
    channelGainOverNoise = channelGaindB - noiseVariancedBm;
    for m = 1:length(Mrange)
        disp([num2str(m) ' antennas out of ' num2str(length(Mrange))]);
        [Hhat,C,tau_p,Rscaled,H] = functionChannelEstimates(R(1:Mrange(m),1:Mrange(m),:,:,:),channelGainOverNoise,nbrOfRealizations,Mrange(m),K,L,p_UL,f);

        [SE_hardening_MR,SE_hardening_RZF,SE_hardening_MMMSE] = functionComputeSE_DL_hardening(H,Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,Mrange(m),K,L,p_UL,p_DL);

        [SE_MR,SE_RZF,SE_MMMSE] = functionComputeSE_DL_estimation(H,Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,Mrange(m),K,L,p_UL,p_DL);

        SE_MR(SE_hardening_MR>SE_MR) = SE_hardening_MR(SE_hardening_MR>SE_MR);
        SE_RZF(SE_hardening_RZF>SE_RZF) = SE_hardening_RZF(SE_hardening_RZF>SE_RZF);
        SE_MMMSE(SE_hardening_MMMSE>SE_MMMSE) = SE_hardening_MMMSE(SE_hardening_MMMSE>SE_MMMSE);

        sumSE_MR(m,1,n) = mean(sum(SE_MR,1));
        sumSE_RZF(m,1,n) = mean(sum(SE_RZF,1));
        sumSE_MMMSE(m,1,n) = mean(sum(SE_MMMSE,1));

            clear H Hhat C Rscaled;
    end
    clear R;
    
end


%% Plot the simulation results
    figure(1);
    hold on; box on;
    
    plot(Mrange,mean(sumSE_MMMSE(:,1,:),3),'rd-','LineWidth',1);
    plot(Mrange,mean(sumSE_RZF(:,1,:),3),'k-.','LineWidth',1);
    plot(Mrange,mean(sumSE_MR(:,1,:),3),'bs-','LineWidth',1);
    
    xlabel('Number of antennas');
    ylabel('Average SE [bit/s/Hz/cell]');
    legend('MMSE','RZF','MR','Location','NorthWest');
    ylim([0 60]);
