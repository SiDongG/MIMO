%% Interfering Power of MS under LoS Massive MIMO
Antenna_size=[1,10,100,1000];

MS_angle=20*pi/180; % Target MS azimuth plane angle
MS_interfere=linspace(-pi,pi); 

d=0.5;

Beta1=1;  %Intracell gain
Beta2=1;  %Intercell gain
Beta=1;   %Large Scale Fading Constant 

g=zeros(length(MS_interfere),length(Antenna_size));

%Channel0=sqrt(Beta)*exp(2*pi*1i*d*(N-1)*sin(MS_angle));  %Target Channel Response
%Channel1=sqrt(Beta)*exp(2*pi*1i*d*(N-1)*sin(theta));     %Interfere Channel Response

for count=1:4
    N=Antenna_size(count);
    j=1;
    Channel0=sqrt(Beta)*exp(2*pi*1i*d*(N-1)*sin(MS_angle)*(0:N-1)');
    while j<101
        theta=MS_interfere(j);
        Channel1=sqrt(Beta)*exp(2*pi*1i*d*(N-1)*sin(theta)*(0:N-1)');
        g(count,j)=abs(Channel0'*Channel1)^2/N;
        j=j+1;
    end
end

figure()
hold on;
plot(MS_interfere,g(1,:),'k-','LineWidth',1);
plot(MS_interfere,g(2,:),'r--','LineWidth',1);
plot(MS_interfere,g(3,:),'b-.','LineWidth',1);
plot(MS_interfere,g(4,:),'g-.','LineWidth',1);

set(gca,'Yscale','log');
xlim([-pi pi]);
ylim([1e-6 1e3]);

legend('M=1','M=10','M=100','M=1000','Location','NorthEast');

xlabel('Azimuth angle of Interfereing MS');
ylabel('Interference Power');