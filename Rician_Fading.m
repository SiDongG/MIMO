%% Calculate the Outage Power under Rician Fading-Sidong Guo

Best_channel_m=1.93*10^-6;
Best_channel_v=1.74*10^-12;
Worst_channel_m=2.41*10^-8;
Worst_channel_v=4.69*10^-16;

Best=zeros(1,26);
Worst=zeros(1,26);  %Pre-allocating for speed

for i=-10:1:15
    i_linear=10^(i/10); %Convert to Linear Scale
    cutoff_best=i_linear*Best_channel_m;
    cutoff_worst=i_linear*Worst_channel_m;
    Shadow_margin_B=Best_channel_m-cutoff_best;
    Shadow_margin_W=Worst_channel_m-cutoff_worst;
    Best(i+11)=qfunc(Shadow_margin_B/sqrt(Best_channel_v));
    Worst(i+11)=qfunc(Shadow_margin_W/sqrt(Worst_channel_v));
end
x=-10:1:15;

figure()
semilogy(x,Best);
hold on 
semilogy(x,Worst);