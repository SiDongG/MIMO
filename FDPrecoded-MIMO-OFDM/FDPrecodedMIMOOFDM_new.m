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
Block_Num=1; %Number of Blocks
SNR=100;
Var_n=1/SNR; %Noise Variance
Ps=10;%Total Power Constraint
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
%% Mapping 
Bits=randi(0:1,[1,Rx*N*Block_Num*log2(M)]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==4
    Bits3=qammod(Bits2,M)*sqrt(0.5);
end
%% Carrier-Wise FFT
Symbol=reshape(Bits3,Tx*N,1,Block_Num);
FFT=dftmtx(N)/sqrt(N);
IFFT=conj(FFT);
Symbol1=zeros(size(Symbol));
for count=1:Block_Num
    Symbol1(:,:,count)=[FFT*Symbol(1:N,:,count);FFT*Symbol(N+1:2*N,:,count)];
end
%% Channel Precoding and Precoding Matrix
%Generate MIMO channel matrix, which is a concatenated 2 dimensional matrix
DD=zeros(N*Rx,N*Tx);
Channel=zeros(1,4,Tx,Rx);
H_bar=zeros(N*Tx,N*Rx);
SS=zeros(N*Tx,N*Rx);
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
        H_bar(N*(i-1)+1:N*i,N*(j-1)+1:N*j)=R*H0*T*IFFT;
        H=zeros(N+L-1,N);
        for count=1:N+L-1
            for m=1:N
                if count-m+1>L
                    H(count,m)=0;
                elseif count-m<0
                    H(count,m)=0;
                else
                    H(count,m)=h(count-m+1);
                end
            end
        end
        H(1:L-1,:)=H(1:L-1,:)+H(N+1:N+L-1,:);
        H=H(1:N,:);
        D=FFT*H*IFFT;
        DD(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=D;
        SS(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=inv(D'*D+Var_n/Ps*eye(N))*D';
    end
end
%% Power Scaling Factor
Sum=trace(SS'*SS);
a=sqrt(Sum/(Tx*N*Ps));
TT=SS/a;
%% Apply Precoding Matrix  
Symbol2=zeros(Tx*N,1,Block_Num);
for count=1:Block_Num
    Symbol2(:,:,count)=TT*Symbol1(:,:,count);
end
%% Apply Channel 
Symbol3=zeros(Tx*N,1,Block_Num);
for count=1:Block_Num
    Symbol3(:,:,count)=H_bar*Symbol2(:,:,count);
end
%% Gain Control
Symbol3=Symbol3*a;
%% Demapping 



