%% Spatial Correlation Matrix  
function R=SpatialCorrelation(M,theta,SD,d)
%Antenna Separation Distance
SD_rad=SD*pi/180;
firstRow = zeros(M,1);
for i = 1:M
    d2 = d*(i-1);
    F = @(Delta)exp(1i*2*pi*d2*sin(theta+Delta)).*exp(-Delta.^2/(2*SD_rad^2))/(sqrt(2*pi)*SD_rad);
    firstRow(i) = integral(F,-20*SD_rad,20*SD_rad);
end

R = toeplitz(firstRow);