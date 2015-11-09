%% same as radius_vs_nbp_data but without using the phages in the model 
% comparison

h = 2.5;
Vbp = (1)^2*(0.34)*pi;
rData = [ 54 85 60 45 61 26 60 90 55 55 50 84 60 87  61]/2;
rin = rData - h;
nbpData = [ 39 100 34 12 45 5 38 121 40  42 35 147 44 146  42]*1000;
% pRvN = polyfit(log(rData),log(nbpData),1);
% nbpFit =exp(pRvN(2))*rData.^(pRvN(1)) ;

%% fiting the data to the formula: nbp = c*rData^3
c = exp(mean(log(nbpData)) -3*mean(log(rin)));
fillfit = c*3*Vbp/(4*pi);
nbpFit2 = c*(rin).^3;
rFitPlot = min(rin):0.01:max(rin);
nbpFitPlot = c*(rFitPlot).^3;

%% error
n = length(rin);
x = log(rin);
y = log(nbpData);
yFit = log(nbpFit2);
VarErr = sum( (log(y) -log(yFit)).^2)/(n-2);
%%
errLogC =sqrt(VarErr*(1/(n-2) + mean(x)^2/sum((x - mean(x)).^2)));
errC = errLogC*c;