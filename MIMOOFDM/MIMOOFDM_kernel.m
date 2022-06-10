clear; clc; close all;
Tx=2; %Number of Transmit Antenna
Rx=2; %Number of Receive Antenna 
L=4;  %Channel Length
C=4;  %CP Length
M=4;  %4-QAM
N=64; %Block Size
Block_Num=Tx*2; %Number of Blocks
loop_Num=100;

%%
total=zeros(1,11,2);  %preallocating for Speed, SNR from 0 to 20
ratio=zeros(1,11,2);
for dB=0:4:40
    disp(dB);
    SNR=10^(dB/10);
    for loop=1:loop_Num
        [Bits_MMSE,Bits_ZF,Bits]=MIMOOFDM1(Tx,Rx,L,C,M,N,Block_Num,SNR);
        ratio(1,dB/4+1,1)=sum(Bits~=Bits_ZF)/(length(Bits));
        ratio(1,dB/4+1,2)=sum(Bits~=Bits_MMSE)/(length(Bits));
        total(1,dB/4+1,1)=total(1,dB/4+1,1)+ratio(1,dB/4+1,1);
        total(1,dB/4+1,2)=total(1,dB/4+1,2)+ratio(1,dB/4+1,2);
    end
end

total=total/loop_Num;
figure()
box on; hold on;
plot(0:4:40,total(:,:,1),'bx-');
plot(0:4:40,total(:,:,2),'rx-');
set(gca,'Yscale','log');
ylim([1e-5 1]);
xlabel('SNR(dB)');
ylabel('Ber');
legend('MIMO-OFDM-ZF','MIMO-OFDM-MMSE')