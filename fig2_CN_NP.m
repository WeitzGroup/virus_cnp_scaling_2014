

lw = 3;
ms = 30;
fs  = 20;
redRat = [106 16 1];   % 
hetRat = [69 16 1];
NP_Rf = redRat(2)/redRat(3); 
CN_Rf = redRat(1)/redRat(2);
CP_Rf = redRat(1)/redRat(3);
CN_het = hetRat(1)/hetRat(2);
NP_het = hetRat(2)/hetRat(3);
[c,n,p] = cnp_head_r(10:150);

[c2,n2,p2] = cnp_head_r([ 50  15  10]);
NP_bp = 7.5/2;
CN_bp = 19.5/7.5;

figure
plot(n./p,c./n,'-r','linewidth',lw)
hold on

plot([NP_bp NP_Rf NP_het]  ,[CN_bp CN_Rf CN_het],'.k','markersize',ms)
plot(n2./p2,c2./n2,'*','markersize',10)
ylim([1 8])
xlim([1 22])

xlabel('N:P','interpreter','latex','fontsize',20)
ylabel('C:N','interpreter','latex','fontsize',20)

printpdf('fig2_CN_NP_42')