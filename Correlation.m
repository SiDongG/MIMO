%% Correlation Coefficient of two Channel estimation of random MS
M=10;
d=0.5;   %Antenna Spacing
SD=10;   %Angular Standard Deviation in degree
MS_angle=90*pi/180;

theta=linspace(-pi,+pi);
theta_degree=theta*180/pi;
SNRdB = 10;
SNR = 10.^(SNRdB/10);
SNRdB2 = 0;
SNR2 = 10.^(SNRdB2/10);

Correlation=zeros(length(theta),1);
R1=SpatialCorrelation(M,MS_angle,SD,d);

for i=1:length(theta)
    R2=SpatialCorrelation(M,theta(i),SD,d);
    R1_matrix=R1(1:M,1:M);
    R2_matrix=R2(1:M,1:M);
    trace1=SNR*abs(trace(R1_matrix*((SNR*R1_matrix+SNR2*R2_matrix+eye(M))\R1_matrix)));
    trace2=SNR2*abs(trace(R2_matrix*((SNR*R1_matrix+SNR2*R2_matrix+eye(M))\R2_matrix)));
    trace3=abs(trace(R1_matrix*((SNR*R1_matrix+SNR2*R2_matrix+eye(M))\R2_matrix)));
    Correlation(i,1)=(sqrt(SNR*SNR2)*trace(3))/sqrt(trace1*trace2);
end

figure;
hold on; box on;

plot(theta_degree,Correlation,'LineWidth',1);

xlabel('Angle of interfering UE [degree]');
ylabel('Antenna-averaged correlation coefficient');
legend('M=100');