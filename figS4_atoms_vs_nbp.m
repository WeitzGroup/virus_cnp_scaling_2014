%clear all
viral_data_full
cnp_virus

mark = {'*','+','x','^','o','s','v'};
lw = 2.5;
fs = 14;
ms = 8;

%%


figure
X = 17;                    %#X = 29.7; A4 paper size
Y = 8;                    %#Y = 21;   A4 paper size
xMargin = 0;               %# left/right margins from page borders
yMargin = 0;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
set(gcf, 'Units','centimeters', 'Position',[5 10 xSize ySize])
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf,'DefaultLineMarkerSize',10);
set(gcf,'DefaultAxesLineWidth',1.5);

axW = 5;
axH = 5;
xpos = 1.5;
ypos = 1.5;

axes('units','centimeters','position',[ xpos ypos axW axH])
plot(nbp,thPhage.C.bp,'-b','linewidth',lw)
hold on
plot(nbp,thPhage.C.bp+ errPhage.C.bp,'--b',nbp,thPhage.C.bp - errPhage.C.bp,'--b','linewidth',lw)
for i= 1:length(phage.name)
plot(phage.length(i),phage.C(i),[ 'k' mark{i}],'markersize',ms,'linewidth',lw)
end

xlabel('Genome size (kbp)','interpreter','latex','fontsize',fs)
ylabel('Atoms in viral head','interpreter','latex','fontsize',fs)
title('Carbon','interpreter','latex','fontsize',fs)
set(gca,'xtick',[10^4:5*10^4:2*10^5],'xticklabel',[10:50:200])
box(gca,'off')
hold off

axes('units','centimeters','position',[ (xpos+axW) ypos axW axH])
plot(nbp,thPhage.N.bp,'-g','linewidth',lw)
hold on
plot(nbp,thPhage.N.bp+ errPhage.N.bp,'--g',nbp,thPhage.N.bp - errPhage.N.bp,'--g','linewidth',lw)
for i= 1:length(phage.name)
plot(phage.length(i),phage.N(i),[ 'k' mark{i}],'markersize',ms,'linewidth',lw)
end

xlabel('Genome size (kbp)','interpreter','latex','fontsize',fs)
title('Nitrogen','interpreter','latex','fontsize',fs)
set(gca,'xtick',[10^4:5*10^4:2*10^5],'xticklabel',[10:50:200])
box(gca,'off')
hold off



axes('units','centimeters','position',[ (xpos+2*axW) ypos axW axH])
plot(nbp,thPhage.P.bp,'-r','linewidth',lw)
hold on
for i= 1:length(phage.name)
plot(phage.length(i),phage.P(i),[ 'k' mark{i}],'markersize',ms,'linewidth',lw)
end

xlabel('Genome size (kbp)','interpreter','latex','fontsize',fs)
title('Phosphorus','interpreter','latex','fontsize',fs)
set(gca,'xtick',[10^4:5*10^4:2*10^5],'xticklabel',[10:50:200])
box(gca,'off')
hold off

%printpdf('figS4_atoms_vs_nbp')
