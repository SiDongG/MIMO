%% Frequency-Domain Precoded MIMO OFDM
clear; clc; close all;
%% Parameter Initialization
% Standard: MMSE Precoding, Complex Orthogonal STBC, 4QAM, TCSI
Tx=2; %Number of Transmit Antenna
Rx=2; %Number of Receive Antenna 
L=4;  %Channel Length
C=4;  %CP Length
M=4;  %4-QAM
N=4; %Block Size
P=N+C;
Block_Num=Rx; %Number of Blocks
SNR=100;
Var_n=1/SNR; %Noise Variance
Ps=10;%Total Power Constraint
%% Mapping 
Bits=randi(0:1,[1,N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Bits3=qammod(Bits2,M)*sqrt(0.5);
end
%% Carrier-wise FFT
Symbol=reshape(Bits3,N,1,Block_Num);
FFT=dftmtx(N)/sqrt(N);
IFFT=conj(FFT);
Symbol1=zeros(size(Symbol));
for count=1:Block_Num
    Symbol1(:,:,count)=FFT*Symbol(:,:,count);
end
%% Channel Precoding
%Generate MIMO channel matrix, which is a concatenated 2 dimensional matrix
DD=zeros(N*Rx,N*Tx);
Channel=zeros(1,4,Tx,Rx);
Ch=zeros(P,P,Tx,Rx);
for i=1:Tx
    for j=1:Rx
        h=(1/sqrt(2*L))*(randn(1,L)+1i*randn(1,L));
        Channel(:,:,i,j)=h;
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
        Ch(:,:,i,j)=H0;
        H=zeros(N+L-1,N);
        for k=1:N+L-1
            for m=1:N
                if k-m+1>L
                    H(k,m)=0;
                elseif k-m<0
                    H(k,m)=0;
                else
                    H(k,m)=h(k-m+1);
                end
            end
        end
        H(1:L-1,:)=H(1:L-1,:)+H(N+1:N+L-1,:);
        H=H(1:N,:);
        DD(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=FFT*H*IFFT;
    end
end
%% Precoding Matrix 
SS=zeros(N*Tx,N*Rx);
SS=inv(DD'*DD+Var_n/Ps*eye(N*Tx))*DD';
%% Power Scaling Factor
sum=0;
for i=1:Tx*N
    for j=1:Rx*N
        sum=sum+abs(SS(i,j))^2;
    end
end

Sum=trace(SS'*SS);
a=sqrt(Sum/(Tx*N*Ps));
TT=SS/a;
%% Apply PreCoding Matrix
Symbol2=zeros(Tx*N,1,Tx,Block_Num);
for k=1:Block_Num
    for i=1:Tx
        for j=1:Rx
            Symbol2(:,:,i,k)=Symbol2(:,:,i,k)+TT(:,N*(j-1)+1:N*j)*Symbol1(:,:,k);
        end
    end
end
%% IFFT 
IFFT2=zeros(N,Tx*N);
for i=1:Tx
    IFFT2(1:N,N*(i-1)+1:N*i)=IFFT;
end
Symbol3=zeros(N,1,Tx,Block_Num);
for count=1:Block_Num
    for i=1:Tx
        Symbol3(:,:,i,count)=IFFT2*Symbol2(:,:,i,count);
    end
end
%% Cyclic Prefix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
Symbol4=zeros(P,1,Tx,Block_Num);
for count=1:Block_Num
    for i=1:Tx
        Symbol4(:,:,i,count)=T*Symbol3(:,:,i,count);
    end
end
%% Channel 
% Not even including the interference channel block because it goes to zero
% anyway
Symbol5=zeros(P,1,Rx,Block_Num);
for count=1:Block_Num
    for i=1:Rx
        for j=1:Tx
            Symbol5(:,:,i,count)=Symbol5(:,:,i,count)+Ch(:,:,j,i)*Symbol4(:,:,j,count);
        end
    end
end
%% Thermal Noise 
Symbol6=zeros(P,1,Rx,Block_Num);
for count=1:Block_Num
    for i=1:Rx       
        nr=randn(P,1);
        ni=randn(P,1);
        Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
        Symbol6(:,:,i,count)=Symbol5(:,:,i,count)+Noise;
    end
end
%% Guard Removal
R=[zeros(N,P-N),eye(N)];
Symbol7=zeros(N,1,Rx,Block_Num);
for count=1:Block_Num
    for i=1:Rx
        Symbol7(:,:,i,count)=R*Symbol6(:,:,i,count);
    end
end
%% Gain Control
Symbol7=Symbol7*a;
%% Demapping 



















