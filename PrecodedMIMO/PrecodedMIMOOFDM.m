%% Precoded MIMO OFDM
function [Bitsre,Bits]=PrecodedMIMOOFDM(Tx,Rx,L,C,M,N,Block_Num,SNR,Eq)
%% Parameter Initialization
% Standard: Complex Orthogonal STBC, 4QAM, CSIR, ZF, MMSE Equalization
% Tx=2; %Number of Transmit Antenna
% Rx=2; %Number of Receive Antenna 
% L=2;  %Channel Length
% C=2;  %CP Length
% M=16;  %16-QAM
% N=4; %Block Size
% Block_Num=Tx*2; %Number of Blocks
% SNR=100;
% Eq=1;
P=N+C;
Var_n=1/SNR; %Noise Variance

%% Mapping 
Bits=randi(0:1,[1,N*Block_Num*log2(M)*2]);
Bits2=zeros(1,length(Bits)/log2(M));
for i=1:log2(M):length(Bits)
    Bits2(1+(i-1)/log2(M))=bin2dec(num2str(Bits(i:i+log2(M)-1)));
end
if M==16
    Symbols=qammod(Bits2,M)*sqrt(1/10);
end
%% Reconfiguring Matrix
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
%% Precoding and AGC
beta=sqrt(Tx/trace(inv(Ch)*(inv(Ch)')));
if Eq==1
    G=beta*inv(Ch);
else
    G=beta*Ch'*inv(Ch*Ch'+Var_n*eye(size(Ch)));
end
Symbols2=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols2(:,:,count)=G*Symbols1(:,:,count);
end
%% Channel
nr=randn(Tx*N,1,Block_Num);
ni=randn(Tx*N,1,Block_Num);
Noise=(1/sqrt(SNR))*(sqrt(2)/2)*(nr+1i*ni);
Symbols3=zeros(N*Tx,1,Block_Num);
for count=1:Block_Num
    Symbols3(:,:,count)=Ch*Symbols2(:,:,count)+Noise(:,:,count);
end
Symbols3=1/beta*Symbols3;
%% Demod 
Bitsre=zeros(size(Bits));
Symbols4=qamdemod(Symbols3/sqrt(1/10),M);
start=1;
for count=1:Block_Num
    for k=1:Tx*N
        dec=dec2bin(Symbols4(k,1,count),log2(M));
        for n=1:length(dec)
            Bitsre(start)=str2double(dec(n));
            start=start+1;
        end
    end
end

