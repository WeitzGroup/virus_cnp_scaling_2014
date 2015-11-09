%%
h = 2.5;

Vbp = (1)^2*(0.34)*pi;
[~,~,r_length_data] = xlsread('phage_genome_radius_with_lipid.xls');
phage_name = r_length_data(2:end,1);
for ii=1:length(phage_name)
    if isa(phage_name{ii},'double')
        phage_name{ii} = int2str(phage_name{ii});
    end
end
rData = cell2mat(r_length_data(2:end,5))./2; 
rin = rData - h;
nbpData = cell2mat(r_length_data(2:end,3));

%%
phage = '44AHJD';
for ii =1:length(phage_name)
    if strcmp(phage_name{ii},phage)
        idx = ii;
        break
    end
end

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
errFillFit = errC*fillfit/c;

%% Residuals
residual = abs(y - yFit);
[~,idxRes] = sort(residual,'descend');

%%
% Plots
figure
%loglog(rData,nbpData,'*',rData,nbpFit2,'-k',rData,nbpFit2,'-r','linewidth',2);
fs = 20;
loglog(rin,nbpData,'ok',rFitPlot,nbpFitPlot,'-k','linewidth',2.5,...
    'markersize',10,'markerfacecolor',[0.1 0.3 0.5]);
xlabel('Capsid internal radius (nm)','fontsize',fs,'interpreter','latex')
ylabel('$n_{bp}$ (number of base pairs)','fontsize',fs,'interpreter','latex')
% hold on

