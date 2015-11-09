clear all
avo = 6.02214179e23;

%% contour  of viral diameter vs viral density.

diam = 30:80/100:100;
vDen = linspace(10^8,10^11,101);
[diam,vDen] = meshgrid(diam,vDen);
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
[~,~,vp] = cnp_head_r(diam/2); 
pDen = vp./avo*10^9.*vDen;
contour(diam,vDen,pDen,[0.01,0.1:0.2:0.7,1,2,3,5, 10:10:40,],'Showtext','on','linewidth',2.5)
set(gca,'yscale','log')
colorbar
%title('Virus P density nmol/L','fontsize',fs)
xlabel('Virus capsid diameter (nm)','fontsize',fs)
ylabel('Virus density (viruses/L)','fontsize',fs) 
%text(5,1.5*10^11,'a','fontsize',25,'fontname','sans')
text(20,1.7*10^11,'a','fontsize',25,'fontname','sans')
printpdf('fig4a_contour_P')
%printpdf('fig4a_contour_P_nocolorbar')