%% 2 by 2 STBC MIMO OFDM
function [Bits_MMSE,Bits_ZF,Bits]=MIMOOFDM1(Tx,Rx,L,C,M,N,Block_Num,SNR)
%% Parameter Initialization
% Standard: Complex Orthogonal STBC, 4QAM, CSIR, ZF, MMSE Equalization
% Tx=2; %Number of Transmit Antenna
% Rx=2; %Number of Receive Antenna 
% L=2;  %Channel Length
% C=2;  %CP Length
% M=4;  %4-QAM
% N=4; %Block Size
% Block_Num=Tx*2; %Number of Blocks
% SNR=100;
P=N+C;
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
FFT=dftmtx(N)/sqrt(N);
IFFT=conj(FFT);
Symbols1=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols1(:,:,count)=[Symbols(:,:,1,count);Symbols(:,:,2,count)];
end
%% Equivalent CFR Channel Matrix 
S=eye(N);
T=[S(2*N-P+1:N,:);S];
R=[zeros(N,P-N),eye(N)];
Ch=zeros(N*Tx,N*Rx);
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
        Ch(N*(j-1)+1:N*j,N*(i-1)+1:N*i)=FFT*R*H0*T*IFFT;
    end
end
%% Channel
nr=randn(Tx*N,1,Block_Num);
ni=randn(Tx*N,1,Block_Num);
Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
Symbols2=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols2(:,:,count)=Ch*Symbols1(:,:,count)+Noise(:,:,count);
end
%% Equalization
G=inv(Ch'*Ch)*Ch';
Symbols_ZF=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols_ZF(:,:,count)=G*Symbols2(:,:,count);
end
G1=inv(Ch'*Ch+Var_n*eye(size(Ch)))*Ch';
Symbols_MMSE=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols_MMSE(:,:,count)=G1*Symbols2(:,:,count);
end
%% Demod 
Bits_ZF=zeros(size(Bits));
Bits_MMSE=zeros(size(Bits));
ZF1=qamdemod(Symbols_ZF/sqrt(1/2),M);
MMSE1=qamdemod(Symbols_MMSE/sqrt(1/2),M);
start=1;
for count=1:Block_Num
    for k=1:Tx*N
        dec=dec2bin(ZF1(k,1,count),log2(M));
        dec2=dec2bin(MMSE1(k,1,count),log2(M));
        for n=1:length(dec)
            Bits_ZF(start)=str2double(dec(n));
            Bits_MMSE(start)=str2double(dec2(n));
            start=start+1;
        end
    end
end




