%% Pilot Contamination with antenna size=100
M=100;
d=0.5;   %Antenna Spacing
SD=10;   %Angular Standard Deviation in degree
MS_angle=20*pi/180;

theta=linspace(-pi,+pi);
theta_degree=theta*180/pi;
SNRdB = [0,20];
SNR1=  10.^(SNRdB/10);
SNRdB2 = 0;
SNR2 = 10.^(SNRdB2/10);

MMSE_correlated = zeros(length(theta),1,length(SNRdB));
MMSE_uncorrelated = zeros(length(theta),1,length(SNRdB));

R1=SpatialCorrelation(M,MS_angle,SD,d);

for i=1:length(theta)
    R2=SpatialCorrelation(M,theta(i),SD,d);
    for j=1:length(SNRdB)
        MMSE_correlated(i,j)=1-SNR2*abs(trace(R1*((SNR2*R1+SNR1(j)*R2+eye(M))\R1)))/trace(R1);
        MMSE_uncorrelated(i,j)=1 - SNR2/(SNR2+SNR1(j)+1);
    end
end

figure;
hold on; box on;

plot(theta_degree,MMSE_correlated(:,1),'k-','LineWidth',1);
plot(theta_degree,MMSE_correlated(:,2),'r-','LineWidth',1);
plot(theta_degree,MMSE_uncorrelated(:,1),'b-','LineWidth',1);
plot(theta_degree,MMSE_uncorrelated(:,2),'y-','LineWidth',1);

xlabel('Interfering Angle');
ylabel('MMSE estimation');


legend('SNRdiff:0dB','SNRdiff:20dB Stronger','SNRdiff:0dB uncorrelated','SNRdiff:20dB uncorrelated');
