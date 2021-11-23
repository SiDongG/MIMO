function [MR,RZF,MMSE] = functionComputeSE_DL_hardening(H,Hhat,C,R,Len_Coherence,Len_pilot,Num_Realizations,M,Num_MS,Num_BS,p_UL,p_DL)

eyeK = eye(Num_MS);
eyeM = eye(M);
C_totM = reshape(p_UL*sum(sum(C,3),4),[M M Num_BS]);

prelogFactor = (Len_Coherence-Len_pilot)./(Len_Coherence);

signal_MR = zeros(Num_MS,Num_BS);
signal_RZF = zeros(Num_MS,Num_BS);
signal_MMMSE = zeros(Num_MS,Num_BS);

interf_MR = zeros(Num_MS,Num_BS);
interf_RZF = zeros(Num_MS,Num_BS);
interf_MMMSE = zeros(Num_MS,Num_BS);


%% Go through all channel realizations
for n = 1:Num_Realizations
    for j = 1:Num_BS
        Hallj = reshape(H(:,n,:,:,j),[M Num_MS*Num_BS]);
        Hhatallj = reshape(Hhat(:,n,:,:,j),[M Num_MS*Num_BS]);
        V_MR = Hhatallj(:,Num_MS*(j-1)+1:Num_MS*j);
        V_RZF = (p_UL*V_MR)/(p_UL*(V_MR'*V_MR)+eyeK);
        V_MMMSE = (p_UL*(Hhatallj*Hhatallj')+C_totM(:,:,j)+eyeM)\(p_UL*V_MR);
        for k = 1:Num_MS            
            if norm(V_MR(:,k))>0
                w = V_MR(:,k)/norm(V_MR(:,k));
                signal_MR(k,j) = signal_MR(k,j) + (w'*H(:,n,k,j,j))/Num_Realizations;
                interf_MR = interf_MR + p_DL*reshape(abs(w'*Hallj).^2,[Num_MS Num_BS])/Num_Realizations;
                w = V_RZF(:,k)/norm(V_RZF(:,k));
                signal_RZF(k,j) = signal_RZF(k,j) + (w'*H(:,n,k,j,j))/Num_Realizations;
                interf_RZF = interf_RZF + p_DL*reshape(abs(w'*Hallj).^2,[Num_MS Num_BS])/Num_Realizations;
                w = V_MMMSE(:,k)/norm(V_MMMSE(:,k));
                signal_MMMSE(k,j) = signal_MMMSE(k,j) + (w'*H(:,n,k,j,j))/Num_Realizations;
                interf_MMMSE = interf_MMMSE + p_DL*reshape(abs(w'*Hallj).^2,[Num_MS Num_BS])/Num_Realizations;
            end
            
        end
        
    end
    
end


%% Prepare to compute the SEs with different pilot reuse factors
MR = zeros(Num_MS,Num_BS,length(Len_Coherence));
RZF = zeros(Num_MS,Num_BS,length(Len_Coherence));
MMSE = zeros(Num_MS,Num_BS,length(Len_Coherence));


for n = 1:length(Len_Coherence)
    MR(:,:,n) = prelogFactor(n)*real(log2(1+(p_DL*abs(signal_MR).^2) ./ (interf_MR - p_DL*abs(signal_MR).^2 + 1)));
    RZF(:,:,n) = prelogFactor(n)*real(log2(1+(p_DL*abs(signal_RZF).^2) ./ (interf_RZF - p_DL*abs(signal_RZF).^2 + 1)));
    MMSE(:,:,n) = prelogFactor(n)*real(log2(1+(p_DL*abs(signal_MMMSE).^2) ./ (interf_MMMSE - p_DL*abs(signal_MMMSE).^2 +1)));
end