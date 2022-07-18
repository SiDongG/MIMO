function [R,channelGaindB]=CellSetup(L,K,M,accuracy,ASDdeg)
%% Inputs
% L=Number of BSs, one per cell
% K=Number of User Equipment per cell
% M=Number of antennas on a BS
% accuracy=Compute exact correlation matrices from the local scattering
% model if approx=1. Compute a small-angle approximation of the model if
% approx=2.
% ASDdeg=Angular standard deviation around the nominal angle
%% Outputs
%R             = M x M x K x L x L matrix with spatial correlation matrices
%                for all UEs in the network. R(:,:,k,j,l) is the correlation
%                matrix for the channel between UE k in cell j and the BS
%                in cell l. This matrix is normalized such that trace(R)=M.
%channelGaindB = K x L x L matrix containing the average channel gains in
%                dB of all the channels. The product 
%                R(:,:,k,j,l)*10^(channelGaindB(k,j,l)/10) is the full
%                spatial channel correlation matrix.
%% Model parameters
%Set the length in meters of the total square area
squareLength = 1000; % 0.25km as length for each rectangular cell, with 16 cells we have 1km in total size. 
%Number of BSs per dimension
nbrBSsPerDim = sqrt(L);
%Pathloss exponent
alpha = 3.76;
%Average channel gain in dB at a reference distance of 1 meter. Note that
%-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
constantTerm = -35.3;
%Standard deviation of shadow fading
sigma_sf = 10;
%Minimum distance between BSs and UEs
minDistance = 35;
%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance
%Distance between BSs in vertical/horizontal direction
interBSDistance = squareLength/nbrBSsPerDim; % Literally the same as 0.25km

%Deploy BSs on the grid
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
% Horizontal
%   125   375   625   875   
%   125   375   625   875
%   125   375   625   875
%   125   375   625   875    
% Vertical
%   125   125   125   125
%   375   375   375   375
%   625   625   625   625
%   875   875   875   875
%Compute all nine alternatives of the BS locations when using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%Prepare to put out UEs in the cells
UEpositions = zeros(K,L);
perBS = zeros(L,1);

%Prepare to store normalized spatial correlation matrices
R = zeros(M,M,K,L,L,length(ASDdeg));

%Prepare to store average channel gain numbers (in dB)
channelGaindB = zeros(K,L,L);


%% Go through all the cells
for n = 1:L
    %Put out K UEs in the cell, uniformly at random. The procedure is
    %iterative since UEs that do not satisfy the minimum distance are
    %replaced with new UEs, Until the number of UEs reach the number of K
    %length(posXY) is either 1 or 0. 
    while perBS(n)<K     
        %Put out new UEs
        UEremaining = K-perBS(n);
        posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posXY = posX + 1i*posY;       
        %Keep those that satisfy the minimum distance
        posXY = posXY(abs(posXY)>=minDistance);        
        %Store new UEs
        UEpositions(perBS(n)+1:perBS(n)+length(posXY),n) = posXY + BSpositions(n);
        perBS(n) = perBS(n)+length(posXY);       
    end  
    %Go through all BSs
    for j = 1:L
        %Compute the distance from the UEs in cell l to BS j with a wrap
        %around topology, where the shortest distance between a UE and the
        %nine different locations of a BS is considered 
        [distancesBSj,whichpos] = min(abs( repmat(UEpositions(:,n),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        %Compute average channel gain using the large-scale fading model in
        %(2.3), while neglecting the shadow fading 
        channelGaindB(:,n,j) = constantTerm - alpha*10*log10(distancesBSj);        
        %Compute nominal angles between UE k in cell l and BS j, and
        %generate spatial correlation matrices for the channels using the
        %local scattering model
        for k = 1:K           
            angleBSj = angle(UEpositions(k,n)-BSpositionsWrapped(j,whichpos(k)));           
            if accuracy == 1 %Use the exact implementation of the local scattering model              
                for spr = 1:length(ASDdeg)                 
                    R(:,:,k,n,j,spr) = functionRlocalscattering(M,angleBSj,ASDdeg(spr),antennaSpacing);                
                end          
            elseif accuracy == 2 %Use the approximate implementation of the local scattering model        
                for spr = 1:length(ASDdeg)                    
                    R(:,:,k,n,j,spr) = functionRlocalscatteringApprox(M,angleBSj,ASDdeg(spr),antennaSpacing);                    
                end                
            end            
        end        
    end        
    %Go through all UEs in cell l and generate shadow fading realizations
    for k = 1:K       
        %Generate shadow fading realizations
        shadowing = sigma_sf*randn(1,1,L);
        channelGainShadowing = channelGaindB(k,n,:) + shadowing;        
        %Check if another BS has a larger average channel gain to the UE
        %than BS l
        while channelGainShadowing(n) < max(channelGainShadowing)            
            %Generate new shadow fading realizations (until all UEs in cell
            %l has its largest average channel gain from BS l)
            shadowing = sigma_sf*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,n,:) + shadowing;            
        end        
        %Store average channel gains with shadowing fading
        channelGaindB(k,n,:) = channelGainShadowing;        
    end    
end
end