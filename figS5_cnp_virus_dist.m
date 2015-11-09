%nucl = fastaread('DNA/dsDNA.fasta');
%nucl_ss = fastaread('DNA/ssDNA.fasta')
load DNA_sequences

mCbp = [19.5 7.5 2];   % C:N:P molecular formula of 2 nucleotides
avo = 6.02214179e23;
atmMass = [12.0107 14.0067 30.973762]; % [C N P]
fgMass = atmMass/avo*10^(15);

protein_char_all
Vbp = (1)^2*(0.34)*pi;   %volume of a base pair
radius_vs_nbp_data  %script to calculate the fit for fillfit
fillC =fillfit;

r_bp =@(x) ((x*Vbp*3)/(4*pi*fillC)).^(1/3)+ h;
rss_bp =@(x) ((x*Vbp*3)/(4*pi*fillC*2)).^(1/3)+ h;  %divide  Vbp/2 for ssDNA 
capsid_carbon = @(x) 4*pi*pCM(1)*Den*(x.^3 - (x-h).^3)/3;
capsid_nitrogen = @(x) 4*pi*pNM(1)*Den*(x.^3 - (x-h).^3)/3;

dsLen = zeros(1,length(nucl));
ssLen = zeros(1,length(nucl_ss));
for ind= 1 :length(nucl);
    dsLen(ind)  = length(nucl(ind).Sequence);
end
for ind = 1 :length(nucl_ss);
    ssLen(ind) = length(nucl_ss(ind).Sequence);
end

dsDNA_C = mCbp(1)*dsLen;
dsDNA_N= mCbp(2)*dsLen;
dsHead_P = mCbp(3)*dsLen;
dsCapsid_C =capsid_carbon(r_bp(dsLen)); 
dsCapsid_N = capsid_nitrogen(r_bp(dsLen));
dsHead_C = dsDNA_C + dsCapsid_C;
dsHead_N = dsDNA_N + dsCapsid_N;

ssDNA_C = mCbp(1)/2*ssLen; %divide mCbp/2 for ssDNA
ssDNA_N= mCbp(2)/2*ssLen;
ssHead_P = mCbp(3)/2*ssLen;
ssCapsid_C =capsid_carbon(rss_bp(ssLen)); 
ssCapsid_N = capsid_nitrogen(rss_bp(ssLen));
ssHead_C = ssDNA_C + ssCapsid_C;
ssHead_N = ssDNA_N + ssCapsid_N;
%% figure including only dsDNA
figure
fs = 14;
sepx = 1.4;
sepy = 1.3;
X = 1.1*8.7;                     %#X = 29.7; A4 paper size
Y = 1.7*X;                  %#Y = 21;   A4 paper size
xMargin = 0;               %# left/right margins from page borders
yMargin = 0;               %# bottom/top margins from page borders
xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
set(gcf, 'Units','centimeters', 'Position',[5 6 xSize ySize])
set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])

ms = 4;
colds =[0.1 0.3 0.5];
colss =[0.5 0.1 0.1];
xl = [0.0002 0.2];
yl = [0 300];

dsC = dsHead_C*fgMass(1);
nBins= ceil(sqrt(length(dsC)));
xC = logspace(log10(min(dsC)),log10(max(dsC)),nBins);
dsHistC =  hist(dsC,xC);
axes('units','centimeters','position',[sepx (sepy+2*(Y-2*sepy)/3)  X-1.5*sepx (Y-2*sepy)/3])
area(xC,dsHistC,-1,'linewidth',2,'linestyle','none','facecolor',colds)
hold on
plot(xC,dsHistC,'o-k','linewidth',2,'markersize',ms,'markerfacecolor',colds)
hold off
set(gca,'xtick',[],'ytick',[0:50:250])
ylim(yl);
set(gca,'xscale','log')
text( 5e-4, 250,'Carbon','interpreter','latex','fontsize',13)
xlim(xl)

dsN = dsHead_N*fgMass(2);
nBins= ceil(sqrt(length(dsN)));
xN = logspace(log10(min(dsN)),log10(max(dsN)),nBins);
dsHistN =  hist(dsN,xN);
axes('units','centimeters','position',[sepx  (sepy+(Y-2*sepy)/3)  (X-1.5*sepx) (Y-2*sepy)/3])
area(xN,dsHistN,-1,'linewidth',2,'linestyle','none','facecolor',colds)
hold on
plot(xN,dsHistN,'o-k','linewidth',2,'markersize',ms,'markerfacecolor',colds)%,xN,ssHistN','ro')
hold off
ylabel('Number of Phages','interpreter','latex','fontsize',fs)
ylim(yl)
set(gca,'xscale','log','xtick',[],'ytick',[0:50:250])
text( 5e-4, 250,'Nitrogen','interpreter','latex','fontsize',13)
xlim(xl)

dsP = dsHead_P*fgMass(2);
nBins= ceil(sqrt(length(dsP)));
xP = logspace(log10(min(dsP)),log10(max(dsP)),nBins);
dsHistP =  hist(dsP,xP);
axes('units','centimeters','position',[sepx sepy (X-1.5*sepx) (Y-2*sepy)/3])
area(xP,dsHistP,-1,'linewidth',2,'linestyle','none','facecolor',colds)
hold on
plot(xP,dsHistP,'o-k','linewidth',2,'markersize',ms,'markerfacecolor',colds)
hold off
ylim(yl)
set(gca,'xscale','log','ytick',[0:50:250])
xlabel('Element in Head (fg)','interpreter','latex','fontsize',fs)
xlim(xl)
text( 5e-4, 250,'Phosphorus','interpreter','latex','fontsize',13)

%printpdf('figS5_cnp_virus_dist')
