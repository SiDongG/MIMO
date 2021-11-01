%% MMSE estimation in terms of receiving SNR 
M=[1,5,10,50,100];
d=0.5;   %Antenna Spacing
SD=10;   %Angular Standard Deviation in degree

theta=linspace(-pi,+pi);
SNRdB = -10:1:20;
SNR = 10.^(SNRdB/10);

MMSE=zeros(length(SNR),length(theta),length(M));


for i=1:length(theta)
    R=SpatialCorrelation(max(M),theta(i),SD,d);
    for j=1:length(SNR)
        for k=1:length(M)
            Size=M(k);
            R_matrix=R(1:Size,1:Size);  %Preallocating R matrix 
            MMSE(j,i,k)=real(trace(R_matrix-SNR(j)*R_matrix*((SNR(j)*R_matrix+eye(Size))\R_matrix)))/trace(R_matrix);
        end
    end
end

figure;
hold on; box on;

plot(SNRdB,mean(MMSE(:,:,1),2),'k-','LineWidth',1);
plot(SNRdB,mean(MMSE(:,:,2),2),'r-','LineWidth',1);
plot(SNRdB,mean(MMSE(:,:,3),2),'b-','LineWidth',1);
plot(SNRdB,mean(MMSE(:,:,4),2),'g-','LineWidth',1);
plot(SNRdB,mean(MMSE(:,:,5),2),'y-','LineWidth',1);



xlabel('Receiving SNR [dB]');
ylabel('MMSE');
set(gca,'YScale','log');

legend('M=1','M=5','M=10','M=50','M=100','Location','SouthWest');