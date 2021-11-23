function [MR,RZF,MMSE] = functionComputeSE_DL_estimation(H_channel,H_channelest,C,~,Len_coherence,Len_pilot,Num_Realizations,M,Num_MS,Num_BS,p_UL,p_DL)

eyeK = eye(Num_MS);
eyeM = eye(M);
Est_err_corr_Sum = reshape(p_UL*sum(sum(C,3),4),[M M Num_BS]);
prelogFactor = (Len_coherence-Len_pilot)./Len_coherence;
gain_MR = zeros(Num_MS,Num_BS,Num_MS);
gain2_MR = zeros(Num_MS,Num_BS,Num_MS);
signal_MR = zeros(Num_MS,Num_BS,Num_Realizations);
intraInterf_MR = zeros(Num_MS,Num_BS,Num_Realizations);
interInterf_MR = zeros(Num_MS,Num_BS);

%% Regularized Zero-forcing Precoding
gain_RZF = zeros(Num_MS,Num_BS,Num_MS);
gain2_RZF = zeros(Num_MS,Num_BS,Num_MS);
signal_RZF = zeros(Num_MS,Num_BS,Num_Realizations);
intraInterf_RZF = zeros(Num_MS,Num_BS,Num_Realizations);
interInterf_RZF = zeros(Num_MS,Num_BS);
%% MMSE Precoding
gain_MMSE = zeros(Num_MS,Num_BS,Num_MS);
gain2_MMSE = zeros(Num_MS,Num_BS,Num_MS);
signal_MMSE = zeros(Num_MS,Num_BS,Num_Realizations);
intraInterf_MMMSE = zeros(Num_MS,Num_BS,Num_Realizations);
interInterf_MMMSE = zeros(Num_MS,Num_BS);


%% Go through all channel realizations
for n = 1:Num_Realizations
    for j = 1:Num_BS
        Hhatallj = reshape(H_channelest(:,n,:,:,j),[M Num_MS*Num_BS]);
        Hallj = reshape(H_channel(:,n,:,:,j),[M Num_MS*Num_BS]);
        V_MR = Hhatallj(:,Num_MS*(j-1)+1:Num_MS*j);
        V_RZF = (p_UL*V_MR)/(p_UL*(V_MR'*V_MR)+eyeK);
        V_MMMSE = (p_UL*(Hhatallj*Hhatallj')+Est_err_corr_Sum(:,:,j)+eyeM)\(p_UL*V_MR);
        for k = 1:Num_MS
            if norm(V_MR(:,k))>0
                w = V_MR(:,k)/norm(V_MR(:,k));
                signal_MR(k,j,n) = abs(w'*H_channel(:,n,k,j,j))^2;
                innerproductsAll = reshape(w'*Hallj,[Num_MS Num_BS]);
                interferenceIntra = Hallj(:,Num_MS*(j-1)+1:Num_MS*j)'*w;
                interInterf_MR = interInterf_MR + abs(innerproductsAll).^2/Num_Realizations;
                interInterf_MR(:,j) = interInterf_MR(:,j) - abs(interferenceIntra).^2/Num_Realizations;
                intraInterf_MR(:,j,n) = intraInterf_MR(:,j,n) + abs(interferenceIntra).^2;
                intraInterf_MR(k,j,n) = intraInterf_MR(k,j,n) - signal_MR(k,j,n);
                gain_MR(:,j,k) = gain_MR(:,j,k) + interferenceIntra/Num_Realizations;
                gain2_MR(:,j,k) = gain2_MR(:,j,k) + abs(interferenceIntra).^2/Num_Realizations;
                
                w = V_RZF(:,k)/norm(V_RZF(:,k));
                signal_RZF(k,j,n) = abs(w'*H_channel(:,n,k,j,j))^2;
                innerproductsAll = reshape(w'*Hallj,[Num_MS Num_BS]);
                interferenceIntra = Hallj(:,Num_MS*(j-1)+1:Num_MS*j)'*w;
                interInterf_RZF = interInterf_RZF + abs(innerproductsAll).^2/Num_Realizations;
                interInterf_RZF(:,j) = interInterf_RZF(:,j) - abs(interferenceIntra).^2/Num_Realizations;
                intraInterf_RZF(:,j,n) = intraInterf_RZF(:,j,n) + abs(interferenceIntra).^2;
                intraInterf_RZF(k,j,n) = intraInterf_RZF(k,j,n) - signal_RZF(k,j,n);
                gain_RZF(:,j,k) = gain_RZF(:,j,k) + interferenceIntra/Num_Realizations;
                gain2_RZF(:,j,k) = gain2_RZF(:,j,k) + abs(interferenceIntra).^2/Num_Realizations;


                w = V_MMMSE(:,k)/norm(V_MMMSE(:,k));
                signal_MMSE(k,j,n) = abs(w'*H_channel(:,n,k,j,j))^2;
                innerproductsAll = reshape(w'*Hallj,[Num_MS Num_BS]);
                interferenceIntra = Hallj(:,Num_MS*(j-1)+1:Num_MS*j)'*w;
                interInterf_MMMSE = interInterf_MMMSE + abs(innerproductsAll).^2/Num_Realizations;
                interInterf_MMMSE(:,j) = interInterf_MMMSE(:,j) - abs(interferenceIntra).^2/Num_Realizations;
                intraInterf_MMMSE(:,j,n) = intraInterf_MMMSE(:,j,n) + abs(interferenceIntra).^2;
                intraInterf_MMMSE(k,j,n) = intraInterf_MMMSE(k,j,n) - signal_MMSE(k,j,n);
                gain_MMSE(:,j,k) = gain_MMSE(:,j,k) + interferenceIntra/Num_Realizations;
                gain2_MMSE(:,j,k) = gain2_MMSE(:,j,k) + abs(interferenceIntra).^2/Num_Realizations;
            end
            
        end
        
    end
    
end


%% Prepare to compute the SEs with different pilot reuse factors
MR = zeros(Num_MS,Num_BS,length(Len_coherence));
RZF = zeros(Num_MS,Num_BS,length(Len_coherence));
MMSE = zeros(Num_MS,Num_BS,length(Len_coherence));
SE_MR_perfect = mean(log2(1+(p_DL*signal_MR) ./ (p_DL*intraInterf_MR + repmat(p_DL*interInterf_MR,[1 1 Num_Realizations]) + 1)),3);

SE_RZF_perfect = mean(log2(1+(p_DL*signal_RZF) ./ (p_DL*intraInterf_RZF + repmat(p_DL*interInterf_RZF,[1 1 Num_Realizations]) + 1)),3);

SE_MMMSE_perfect = mean(log2(1+(p_DL*signal_MMSE) ./ (p_DL*intraInterf_MMMSE + repmat(p_DL*interInterf_MMMSE,[1 1 Num_Realizations]) + 1)),3);

for n = 1:length(Len_coherence)    
    MR(:,:,n) = prelogFactor(n)*SE_MR_perfect - log2(prod(1+p_DL*(Len_coherence(n)-Len_pilot)*(gain2_MR-abs(gain_MR).^2),3))/Len_coherence(n);
    RZF(:,:,n) = prelogFactor(n)*SE_RZF_perfect - log2(prod(prod(1+p_DL*(Len_coherence(n)-Len_pilot)*(gain2_RZF-abs(gain_RZF).^2),3),4))/Len_coherence(n);    
    MMSE(:,:,n) = prelogFactor(n)*SE_MMMSE_perfect - log2(prod(prod(1+p_DL*(Len_coherence(n)-Len_pilot)*(gain2_MMSE-abs(gain_MMSE).^2),3),4))/Len_coherence(n);
end

MR(MR<0) = 0;
RZF(RZF<0) = 0;
MMSE(MMSE<0) = 0;

