
% main data goes here
% Virus based on this doi:10.1016/j.jmb.2007.02.046 and the Sullivan data
% This data is in femtograms
med4.c=[45.8 60.9];
med4.n=[9.4 9.6];
med4.p=[0.98 0.34]; %[replete limited]
w8012.c=[92.4 132]; 
w8012.n=[20 20.6];
w8012.p=[1.84 0.47];
w8013.c=[213 244];
w8013.n=[50.2 39.8];
w8013.p=[3.34 0.81];
v.bp=50000;
v.r=35; % nm
v.burst=40;
% This data is in molecules
[v.cmol,v.nmol,v.pmol] = cnp_head_r(v.r);
v.c=v.cmol/(6.02*10^23)*12.01*10^15;  % fg
v.n=v.nmol/(6.02*10^23)*14.01*10^15;  % fg
v.p=v.pmol/(6.02*10^23)*30.97*10^15;  % fg

%% plots
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
fs = 11;
lw =3;

axes('units','centimeters','position',[ xpos ypos axW axH])
x=[0:0.001:1];  % From replete to limiting
iso=med4;
lys_med4=lysate_calc(iso,v,x);
tmph=plot(1-x,lys_med4.c./lys_med4.isoc,'b-');
set(tmph,'linewidth',lw);
hold on
tmph=plot(1-x,lys_med4.n./lys_med4.ison,'g-');
set(tmph,'linewidth',lw);
tmph=plot(1-x,lys_med4.p./lys_med4.isop,'r-');
set(tmph,'linewidth',lw);
ylim([0 1]);
set(gca,'fontsize',fs,'xtick',[0 0.2 0.4 0.6 0.8]);
xlabel({'Host conditions ','(0 - limiting, 1 - replete)'},'fontsize',fs,'verticalalignment','top');
ylabel('Fraction of host element in lysate','fontsize',fs,'verticalalignment','bottom');
title({'MED4 w/70nm podo', 'burst size 40'},'fontsize',fs)
% legend
% tmplh = legend('stuff',...);
tmplh = legend('C','N','P',4);
% remove box
% set(tmplh,'visible','off')
legend('boxoff');


axes('units','centimeters','position',[ (xpos+axW) ypos axW axH])
x=[0:0.001:1];  % From replete to limiting
iso=w8012;
lys_w8012=lysate_calc(iso,v,x);
tmph=plot(1-x,lys_w8012.c./lys_w8012.isoc,'b-');
set(tmph,'linewidth',lw);
hold on
tmph=plot(1-x,lys_w8012.n./lys_w8012.ison,'g-');
set(tmph,'linewidth',lw);
tmph=plot(1-x,lys_w8012.p./lys_w8012.isop,'r-');
set(tmph,'linewidth',lw);
ylim([0 1]);
set(gca,'fontsize',fs);
xlabel({'Host conditions ','(0 - limiting, 1 - replete)'},'fontsize',fs,'verticalalignment','top');
%ylabel('Fraction of host element in lysate','fontsize',20,'verticalalignment','bottom');
title({'WH8012 w/70nm podo', 'burst size 40'},'fontsize',fs)
% legend
% tmplh = legend('stuff',...);
tmplh = legend('C','N','P',4);
% remove box
% set(tmplh,'visible','off')
set(gca,'yticklabel',[],'xtick',[0 0.2 0.4 0.6 0.8]);
legend('boxoff');

x=[0:0.001:1];  % From replete to limiting


axes('units','centimeters','position',[ (xpos+2*axW) ypos axW axH])
iso=w8013;
lys_w8013=lysate_calc(iso,v,x);
tmph=plot(1-x,lys_w8013.c./lys_w8013.isoc,'b-');
set(tmph,'linewidth',lw);
hold on
tmph=plot(1-x,lys_w8013.n./lys_w8013.ison,'g-');
set(tmph,'linewidth',lw);
tmph=plot(1-x,lys_w8013.p./lys_w8013.isop,'r-');
set(tmph,'linewidth',lw);
ylim([0 1]);
set(gca,'fontsize',fs);
xlabel({'Host conditions ','(0 - limiting, 1 - replete)'},'fontsize',fs,'verticalalignment','top');
%ylabel('Fraction of host element in lysate','fontsize',20,'verticalalignment','bottom');
title({'WH8013 w/70nm podo', 'burst size 40'},'fontsize',fs)
% Save the data
% save cnp_dist2_analysis
% legend
% tmplh = legend('stuff',...);
tmplh = legend('C','N','P',4);
% remove box
% set(tmplh,'visible','off')
legend('boxoff');
set(gca,'yticklabel',[],'xtick',[0 0.2 0.4 0.6 0.8]);

printpdf('fig3b_stoichiometric_mismatch')
