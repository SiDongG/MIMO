function [R,channelGaindB] = functionExampleSetup(Num_BS,Num_MS,M,accuracy,SD)

%% Model parameters
squareLength = 1000;
nbrBSsPerDim = sqrt(Num_BS);
alpha = 3.76;
constantTerm = -35.3;
sigma_sf = 10;
minDistance = 35;
d = 1/2;
interBSDistance = squareLength/nbrBSsPerDim;
locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[Num_BS 1]);
UEpositions = zeros(Num_MS,Num_BS);
perBS = zeros(Num_BS,1);
R = zeros(M,M,Num_MS,Num_BS,Num_BS,length(SD));
channelGaindB = zeros(Num_MS,Num_BS,Num_BS);


%% Go through all the cells
for l = 1:Num_BS
    while perBS(l)<Num_MS
        UEremaining = Num_MS-perBS(l);
        posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
        posXY = posX + 1i*posY;
        posXY = posXY(abs(posXY)>=minDistance);
        UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
        perBS(l) = perBS(l)+length(posXY);
        
    end
    for j = 1:Num_BS
        [distancesBSj,whichpos] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[Num_MS 1]) ),[],2);
        channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
        for k = 1:Num_MS            
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,whichpos(k)));            
            if accuracy == 1                
                for spr = 1:length(SD)                    
                    R(:,:,k,l,j,spr) = functionRlocalscattering(M,angleBSj,SD(spr),d);                   
                end                
            elseif accuracy == 2                
                for spr = 1:length(SD)                    
                    R(:,:,k,l,j,spr) = functionRlocalscatteringApprox(M,angleBSj,SD(spr),d);                    
                end                
            end            
        end        
    end
    for k = 1:Num_MS
        shadowing = sigma_sf*randn(1,1,Num_BS);
        channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        while channelGainShadowing(l) < max(channelGainShadowing)
            shadowing = sigma_sf*randn(1,1,Num_BS);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            
        end
        channelGaindB(k,l,:) = channelGainShadowing;        
    end    
end