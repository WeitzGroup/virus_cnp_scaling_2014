clear all
radius_vs_nbp_data  %script to calculate the fit for fillfit
fillC =roundn(fillfit,-2); % fraction of the capsid fill with DNA rounded to two decimals
h = 2.5;
nbpR = @(x)  4/3*pi*(x-h).^3./Vbp.*fillC; 

avo = 6.02214179e23;

%% contour  of viral diameter vs viral density.

nbpStart = nbpR(15);
nbpFinal = nbpR(50);
nbp = nbpStart:(nbpFinal-nbpStart)/100:nbpFinal;
vDen = linspace(10^8,10^11,101);
[nbp,vDen] = meshgrid(nbp,vDen);

fs = 15;
figure
X = 16;                    %#X = 29.7; A4 paper size
Y = 12;                    %#Y = 21;   A4 paper size
xMargin = 0;               %# left/right margins from page borders
yMargin = 0;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
set(gcf, 'Units','centimeters', 'Position',[5 5 xSize ySize])
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])

vp = 2*nbp; 
pDen = vp./avo*10^9.*vDen;
contour(nbp,vDen,pDen,[0.01,0.1:0.2:0.7,1,2,3,5, 10:10:40,],'Showtext','on','linewidth',2.5)

set(gca,'yscale','log')
hcb = colorbar;
%title('Virus P density nmol/L','fontsize',fs)
xlabel('Virus genome length (bp)','fontsize',fs)
ylabel('Virus density (viruses/L)','fontsize',fs) 
text(-40000,1.7*10^11,'b','fontsize',25,'fontname','sans')
%printpdf('fig4b_contour_P_nbp')
