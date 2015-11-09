%%gets radius length and cnp composition of the complete head (including
%%decorations) of  the phages in phage_full_collab.xls

mCbp = [19.5 7.5 2]; 
matlabDir = dir('*.txt');
nFiles = length(matlabDir);
fileNames = {};

phage.length = [];
capsid.C = [];
capsid.N =[];
capsid.P =[];
DNA.C = [];
DNA.N =[];
DNA.P =[];

[~,~,rdata] = xlsread('phage_full_collab.xls'); % read the excel file with the info
phage.name = rdata(2:end,1)';       %extract the names
phage.name(cellfun(@(x) any(isnan(x)),phage.name))=[];    %delete [NaN]s

phage.radius =rdata(2:end,6)'; %extract the radiui
phage.radius(cellfun(@(x) any(isnan(x)),phage.radius))=[]; %clean cell
phage.radius = cell2mat(phage.radius); % transform cell to matrix

%% capsid molecular composition
for i= 1:length(phage.name)   
    molCompCap = get_cap_mol_comp(rdata,phage.name{i}); %[C,H,N,O,S]
    capsid.C(i) = molCompCap(1);
    capsid.N(i) = molCompCap(3);
end

%% DNA molecular composition

for i= 1:length(phage.name)  
%     [row,column]= find(strcmp(rdata,phage.name{i}));   %find the row of the phage
%     DNAlength = rdata{row,5};     
%     molCompDNA = DNAlength*mCbp;
    filename = [phage.name{i} '_DNA.txt'];
    nucleobase = fastaread(filename);
    phage.length(i) = length(nucleobase.Sequence);
    molCompDNA = get_DNA_mol_form(phage.name{i});
    DNA.C(i) = molCompDNA(1);
    DNA.N(i) = molCompDNA(2);
    DNA.P(i) = molCompDNA(3);
end

phage.C=capsid.C+DNA.C; phage.N=capsid.N+DNA.N;  phage.P=DNA.P;
        

