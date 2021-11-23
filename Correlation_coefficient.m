%% Correlation Coefficient of two Channel estimation of random MS
M=100;
d=0.5;   %Antenna Spacing
SD=10;   %Angular Standard Deviation in degree
MS_angle=20*pi/180;

theta=linspace(-pi,+pi);
theta_degree=theta*180/pi;
SNRdB = [0,20,40];
SNRdB2 = 0;
SNR2 = 10.^(SNRdB2/10);

Correlation=zeros(length(theta),1,length(SNRdB));
R1=SpatialCorrelation(M,MS_angle,SD,d);

for i=1:length(theta)
    R2=SpatialCorrelation(M,theta(i),SD,d);
    R1_matrix=R1(1:M,1:M);
    R2_matrix=R2(1:M,1:M);
    for j=1:length(SNRdB)
        SNR1=10.^(SNRdB(j)/10);
        Correlation(i,1,j)=sqrt(SNR1*SNR2)*abs(trace(R1_matrix*((SNR1*R1_matrix+SNR2*R2_matrix+eye(M))\R2_matrix)))/sqrt(SNR1*SNR2*abs(trace(R1_matrix*((SNR1*R1_matrix+SNR2*R2_matrix+eye(M))\R1_matrix)))*abs(trace(R2_matrix*((SNR1*R1_matrix+SNR2*R2_matrix+eye(M))\R2_matrix))));
    end
end

figure;
hold on; box on;

plot(theta_degree,Correlation(:,1),'k-','LineWidth',1);
plot(theta_degree,Correlation(:,2),'b-','LineWidth',1);
plot(theta_degree,Correlation(:,3),'r-','LineWidth',1);

xlabel('Intefereing Angle');
ylabel('Correlation Coefficient');
legend('diff(SNR)=0','diff(SNR)=20','diff(SNR)=40');

