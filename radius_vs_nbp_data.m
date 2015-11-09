%%
h = 2.5;
Vbp = (1)^2*(0.34)*pi;
rData = [63 54 85 60 45 61 26 60 90 61 55 55 50 84 70 60 60 87 66 61]/2;
rin = rData - h;
nbpData = [49 39 100 34 12 45 5 38 121 40 40 42 35 147 71 46 44 146 40 42]*1000;
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