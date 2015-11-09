clear all
h = 2.5;
Vbp = (1)^2*(0.34)*pi;   %volume of a base pair
radius_vs_nbp_data

% Plots
figure
loglog(rData,nbpData,'*',rData,nbpFit2,'-k',rData,nbpFit2,'-r','linewidth',2);
fs = 20;
loglog(rin,nbpData,'ok',rFitPlot,nbpFitPlot,'-k','linewidth',2.5,...
    'markersize',10,'markerfacecolor',[0.1 0.3 0.5]);
xlabel('Capsid internal radius (nm)','fontsize',fs,'interpreter','latex')
ylabel('$n_{bp}$ (number of base pairs)','fontsize',fs,'interpreter','latex')

%fitFormula = ['$n_{bp}=' num2str(exp(pRvN(2))) 'r^{' num2str(pRvN(1)) '}$'];
fitFormula2 = ['$n_{bp}=' num2str(c) 'r^3$'];
%legend({'Data from DePaepe 2009',fitFormula2},'interpreter','latex','fontsize',10) 
%printpdf('length_vs_r')