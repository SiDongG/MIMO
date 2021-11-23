function [Hhat_MMSE,C_MMSE,Len_pilot,R,H,Hhat_EW_MMSE,C_EW_MMSE,Hhat_LS,C_LS] = functionChannelEstimates(R,channelGaindB,Num_Realizations,M,Num_MS,Num_BS,p_UL,f)

%% Generate channel realizations
H = (randn(M,Num_Realizations,Num_MS,Num_BS,Num_BS)+1i*randn(M,Num_Realizations,Num_MS,Num_BS,Num_BS));
betas = zeros(Num_MS,Num_BS,Num_BS);
for j = 1:Num_BS
    for l = 1:Num_BS     
        for k = 1:Num_MS       
            if channelGaindB(k,j,l)>-Inf
                betas(k,j,l) = 10^(channelGaindB(k,j,l)/10);
                R(:,:,k,j,l) = betas(k,j,l) * R(:,:,k,j,l);
                Rsqrt = sqrtm(R(:,:,k,j,l));
                H(:,:,k,j,l) = sqrt(0.5)*Rsqrt*H(:,:,k,j,l);
            else                
                betas(k,j,l) = 0;
                R(:,:,k,j,l) = 0;
                H(:,:,k,j,l) = 0;               
            end           
        end       
    end   
end

%% Channel estimation
Len_pilot = f*Num_MS;  
pilotPattern = kron(ones(2,1),[1; 2; 1; 2; 3; 4; 3; 4]);
eyeM = eye(M);
Np = sqrt(0.5)*(randn(M,Num_Realizations,Num_MS,Num_BS,f) + 1i*randn(M,Num_Realizations,Num_MS,Num_BS,f));
Hhat_MMSE = zeros(M,Num_Realizations,Num_MS,Num_BS,Num_BS);
C_MMSE = zeros(M,M,Num_MS,Num_BS,Num_BS);
if nargout >= 5
    Hhat_EW_MMSE = zeros(M,Num_Realizations,Num_MS,Num_BS,Num_BS);
    C_EW_MMSE = zeros(M,M,Num_MS,Num_BS,Num_BS); 
end
if nargout >= 7
    Hhat_LS = zeros(M,Num_Realizations,Num_MS,Num_BS,Num_BS);
    C_LS = zeros(M,M,Num_MS,Num_BS,Num_BS);   
end


%% Go through all cells
for j = 1:Num_BS
    for g = 1:f
        groupMembers = find(g==pilotPattern)';
        yp = sqrt(p_UL)*Len_pilot*sum(H(:,:,:,g==pilotPattern,j),4) + sqrt(Len_pilot)*Np(:,:,:,j,g);
        for k = 1:Num_MS
            PsiInv = (p_UL*Len_pilot*sum(R(:,:,k,g==pilotPattern,j),4) + eyeM);
            if nargout >= 5
                PsiInvDiag = diag(PsiInv);
            end
            for l = groupMembers
                RPsi = R(:,:,k,l,j) / PsiInv;
                Hhat_MMSE(:,:,k,l,j) = sqrt(p_UL)*RPsi*yp(:,:,k);
                C_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - p_UL*Len_pilot*RPsi*R(:,:,k,l,j);
                if nargout >= 5
                    A_EW_MMSE = diag(sqrt(p_UL)*diag(R(:,:,k,l,j)) ./ PsiInvDiag);
                    Hhat_EW_MMSE(:,:,k,l,j) = A_EW_MMSE*yp(:,:,k);
                    productAR = A_EW_MMSE * R(:,:,k,l,j);
                    C_EW_MMSE(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p_UL)*Len_pilot + Len_pilot*A_EW_MMSE*PsiInv*A_EW_MMSE';
                
                end
                if nargout >= 7
                    A_LS = 1/(sqrt(p_UL)*Len_pilot);
                    Hhat_LS(:,:,k,l,j) = A_LS*yp(:,:,k);
                    productAR = A_LS * R(:,:,k,l,j);
                    C_LS(:,:,k,l,j) = R(:,:,k,l,j) - (productAR + productAR') * sqrt(p_UL)*Len_pilot + Len_pilot*A_LS*PsiInv*A_LS';                
                end                
            end            
        end
    end    
end