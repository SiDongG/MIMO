%% 2 by 2 STBC MIMO OFDM
clear; clc; close all;
%% Parameter Initialization
% Standard: Complex Orthogonal STBC, 4QAM, CSIR, ZF, MMSE Equalization
Tx=2; %Number of Transmit Antenna
Rx=2; %Number of Receive Antenna 
L=2;  %Channel Length
C=2;  %CP Length
M=4;  %4-QAM
N=4; %Block Size
P=N+C;
Block_Num=Tx*2; %Number of Blocks
SNR=100;
Var_n=1/SNR; %Noise Variance

%% Mapping 
Bits=randi(0:1,[1,N*Block_Num*log2(M)*2]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Symbols=qammod(Bits2,M)*sqrt(0.5);
end
%% IFFT
Symbols=reshape(Symbols,[N,1,Tx,Block_Num]);
Symbols0=zeros(size(Symbols));
FFT=dftmtx(N)/sqrt(N);
IFFT=conj(FFT);
for count=1:Block_Num
    for i=1:Tx
        Symbols0(:,:,i,count)=FFT*Symbols(:,:,i,count);
    end
end
%% Cyclic Prefix 
Symbols1=zeros(P,1,Tx,Block_Num);
S=eye(N);
T=[S(2*N-P+1:N,:);S];
for count=1:Block_Num
    for i=1:Tx
        Symbols1(:,:,i,count)=T*Symbols0(:,:,i,count);
    end
end
%% Channel and Equalization Matrix 
Ch=zeros(P*Tx,P*Rx);
% G=zeros(P*Tx,P*Rx);
for i=1:Tx
    for j=1:Rx
        h=(1/sqrt(2*L))*(randn(1,L)+1i*randn(1,L));
        a=1;
        H0=zeros(P);
        while a<P+1  %generate the channel matrces
            b=1;
            while b<P+1
                if a-b<0 || a-b>L-1
                    H0(a,b)=0;
                else
                    H0(a,b)=h(a-b+1);
                end
                b=b+1;
            end
            a=a+1;
        end
        Ch(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=H0;
    end
end
% %% Channel Transmission
% Symbols2=zeros(P,1,Rx,Block_Num);
% for count=1:Block_Num
%     for i=1:Rx
%         for j=1:Tx
%             Symbols2(:,:,i,count)=Symbols2(:,:,i,count)+Ch(:,:,j,i)*Symbols1(:,:,j,count);
%         end
%     end
% end
% %% Truncation
% R=[zeros(N,P-N),eye(N)];
% Symbols3=zeros(N,1,Rx,Block_Num);
% for count=1:Block_Num
%     for i=1:Rx
%         Symbols3(:,:,i,count)=R*Symbols2(:,:,i,count);
%     end
% end
% %% FFT
% Symbols4=zeros(N,1,Rx,Block_Num);
% for count=1:Block_Num
%     for i=1:Rx
%         Symbols4(:,:,i,count)=FFT*Symbols3(:,:,i,count);
%     end
% end
% %% Equalization Matrix
% G=zeros(N,N,Rx,Tx);
% for i=1:Rx
%     for j=1:Tx
%         D=FFT*R*Ch(:,:,j,i)*T*IFFT;
%         G(:,:,i,j)=pinv(D);
%     end
% end
% %% Equalization
% Symbols5=zeros(N,1,Rx,Block_Num);
% for count=1:Block_Num
%     for i=1:Rx
%         for j=1:Rx
%             Symbols5(:,:,i,count)=Symbols5(:,:,i,count)+G(:,:,i,j)*Symbols4(:,:,j,count);
%         end
%     end
% end









