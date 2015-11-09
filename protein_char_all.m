avo = 6.02214179e23;
atmWeight = [12.0107 1.00794 14.0067 15.9994 32.065];
h = 2.5;  %(nm) average protein height   
Psv = 0.73e21;
errPsv = 0.03e21;
Den = 1/Psv;
errDen = 1/Psv^2*errPsv;

%cutoffLength =6000;  %getting rid of two outliers.
AA = fastaread('protein\protein_all.fasta');
nProt = length(AA);

protLength = zeros(1,nProt);
protWeight = zeros(1,nProt);
protC = zeros(1,nProt);
protN =  zeros(1,nProt);
molComp = zeros(nProt,5);

    for i = 1:nProt;
    sequence = AA(i).Sequence;
        protLength(i) = length(sequence);
        molComp(i,:) = protein_mol_form(sequence);  %[C H N O S]
        protC(i) = molComp(i,1);
        protN(i) = molComp(i,3);
        protWeight(i) = molComp(i,:)*atmWeight'; 
end
 %%  fitting the data to a line to obtain (# molecules/ unit mass)

 ind = find(protWeight>4e5);
 protWeight(ind) = [];
 protC(ind) = [];
 protN(ind) = [];
 
   pCM = polyfit(protWeight,protC,1);
   pNM = polyfit(protWeight,protN,1);
        
   fCM = polyval(pCM,protWeight);
   fNM = polyval(pNM,protWeight);
    
   pCM(1) = pCM(1)*avo;
   pNM(1) = pNM(1)*avo;
    
  errpCM=sqrt(sum((protC-fCM).^2)/sum((protWeight/avo-mean(protWeight/avo)).^2)...
             /(length(protWeight)-2));
   
  errpNM=sqrt(sum((protN-fNM).^2)/sum((protWeight/avo-mean(protWeight/avo)).^2)...
             /(length(protWeight)-2));
         
dC = pCM(1)*Den;
dN = pNM(1)*Den;

errdC = sqrt((pCM(1)*errDen)^2 + (Den*errpCM)^2);
errdN = sqrt((pNM(1)*errDen)^2 + (Den*errpNM)^2);

%%
% figure
% fs = 15;
% X = 30;                    %#X = 29.7; A4 paper size
% Y = 10;                    %#Y = 21;   A4 paper size
% xMargin = 0;               %# left/right margins from page borders
% yMargin = 0;               %# bottom/top margins from page borders
% xSize = X - 2*xMargin;     %# figure size on paper (widht & hieght)
% ySize = Y - 2*yMargin;     %# figure size on paper (widht & hieght)
% set(gcf, 'Units','centimeters', 'Position',[55 10 xSize ySize])
% set(gcf, 'PaperUnits','centimeters')
% set(gcf, 'PaperSize',[X Y])
% set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
% 
% subplot(1,2,1)
%     plot(protWeight,protC,'o','markersize',10)
%     hold on
%     plot(protWeight,fCM,'-k','linewidth',2)
%     hold off
%     ylabel('Number of Carbon atoms','interpreter','latex','fontsize',fs)
%     xlabel('Mass of protein (u)','interpreter','latex','fontsize',fs)
% 
%     textcorner('a',0.5,0.5,1.5*fs)
% subplot(1,2,2) 
%     plot(protWeight,protN,'o','markersize',10)
%     hold on
%     plot(protWeight,fNM,'-k','linewidth',2)
%     hold off
%     ylabel('Number of Nitrogen atoms','interpreter','latex','fontsize',fs)
%     xlabel('Mass of protein (u)','interpreter','latex','fontsize',fs)
%     textcorner('b',0.5,0.5,1.5*fs)
% 
% printpdf('atomsCN_vs_mass')