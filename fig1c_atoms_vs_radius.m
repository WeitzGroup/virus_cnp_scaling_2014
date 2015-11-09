%clear all
viral_data_full
cnp_virus

%%
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
plot(r,thPhage.C.r,'-b','linewidth',lw)
hold on
plot(r,thPhage.C.r+ errPhage.C.r,'--b',r,thPhage.C.r - errPhage.C.r,'--b','linewidth',lw)
for i= 1:length(phage.name)
plot(phage.radius(i),phage.C(i),[ 'k' mark{i}],'markersize',ms,'linewidth',lw)
end

yl = ylim*0.8;
ylim([ 0 yl(2)])
xlim([0 60])
xlabel('Head radius (nm)','interpreter','latex','fontsize',fs)
ylabel('Atoms in virus head','interpreter','latex','fontsize',fs)
title('Carbon','interpreter','latex','fontsize',fs)
set(gca,'xtick',[0 20 40],'xticklabel',[0 20 40])
box(gca,'off')
hold off

axes('units','centimeters','position',[ (xpos+axW) ypos axW axH])
plot(r,thPhage.N.r,'-g','linewidth',lw)
hold on
%plot(phage.radius,phage.N,'.k','markersize',ms)
plot(r,thPhage.N.r+ errPhage.N.r,'--g',r,thPhage.N.r - errPhage.N.r,'--g','linewidth',lw)
for i= 1:length(phage.name)
plot(phage.radius(i),phage.N(i),[ 'k' mark{i}],'markersize',ms,'linewidth',lw)
end
yl = ylim;
ylim([ 0 yl(2)])
xlim([0 60])
xlabel('Head radius (nm)','interpreter','latex','fontsize',fs)
title('Nitrogen','interpreter','latex','fontsize',fs)
set(gca,'xtick',[20 40],'xticklabel',[ 20 40],'box','off')
ytick = get(gca,'ytick');
ytick = ytick(2:end-1);
set(gca,'ytick',ytick)

hold off

axes('units','centimeters','position',[ (xpos+2*axW) ypos axW axH])
plot(r, thPhage.P.r, '-r', 'linewidth',lw)
hold on
plot(r,thPhage.P.r+ errPhage.P.r,'--r',r,thPhage.P.r - errPhage.P.r,'--r','linewidth',lw)
for i= 1:length(phage.name)
plot(phage.radius(i),phage.P(i),[ 'k' mark{i}],'markersize',ms,'linewidth',lw)
end
yl = ylim;
ylim([ 0 yl(2)])
xlim([0 60])
xlabel('Head radius (nm)','interpreter','latex','fontsize',fs)
title('Phosphorus','interpreter','latex','fontsize',fs)
set(gca,'xtick',[20 40],'xticklabel',[ 20 40],'box','off')
ytick = get(gca,'ytick');
ytick = ytick(2:end-1);
set(gca,'ytick',ytick);
hold off

%printpdf('fig_1c_atoms_vs_radius')
