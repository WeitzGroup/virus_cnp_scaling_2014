load DNA_sequences

%% DNA primary structure
base = 'ATCG';
nucFreqM = zeros(length(nucl),4);
for i =1:length(nucl);
    sequence = nucl(i).Sequence;
    nucFreq = zeros(1,4);
    for ind=1:4
        nucFreq(ind) = length(find(sequence==base(ind)));
    end
    nucFreq =  nucFreq/length(sequence);
    nucFreqM(i,:) = nucFreq;
end

%%
at = sum(nucFreqM(:,[1 2]),2);
cg = sum(nucFreqM(:,[3 4]),2);
A = [5 5 ];
T = [5 2 ];
C = [4 3 ];
G = [5 5 ];
S = [10 0];
m = nucFreqM;
atComp = at*(A+T);
cgComp = cg*(C+G);
molCompNuc = mean(at)*(A+T)+ mean(cg)*(C +G);

a= atComp + cgComp + ones(size(at))*S; 

figure
fs = 14;
X = 30;                     %#X = 29.7; A4 paper size
Y = 7;                  %#Y = 21;   A4 paper size
xMargin = 0;               %# left/right margins from page borders
yMargin = 0;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
set(gcf, 'Units','centimeters', 'Position',[5 8 xSize ySize])
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])

subplot(1,3,2)
hist(a(:,1))
xlabel('Avg. Carbon in Nucleotide','fontsize',fs,'interpreter','latex')
xl = xlim;
yl = ylim;
w = xl(2)-xl(1);
text(xl(1) +0.03*w,yl(2),'\textbf{b}','fontsize',20,'interpreter','latex') 
%title(['mean = ' int2str(mean(a(:,1)))]) 

subplot(1,3,3)
hist(a(:,2))
xlabel('Avg. Nitrogen in Nucleotide','fontsize',fs,'interpreter','latex')
xl = xlim;
yl = ylim;
w = xl(2)-xl(1);
text(xl(1) -0.2*w,yl(2),'\textbf{c}','fontsize',20,'interpreter','latex') 

subplot(1,3,1)
hist(cg)
xlabel('Fraction of C and G','fontsize',fs,'interpreter','latex')
ylabel('Number of sequences','fontsize',fs,'interpreter','latex')
xl = xlim;
yl = ylim;
w = xl(2)-xl(1);
text(xl(1) - 0.2*w,yl(2),'\textbf{a}','fontsize',20,'interpreter','latex') 

%printpdf('hist_nucl_CN')
