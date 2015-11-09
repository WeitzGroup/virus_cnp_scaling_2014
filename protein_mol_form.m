function molCompP = protein_mol_form(seq)
%molComP = [C H N O S]
aminoacids = 'ACDEFGHIKLMNPQRSTVWY';
amiComp = [3 7 1 2 0; 3 7 1 2 1; 4 7 1 4 0;5 9 1 4 0;9 11 1 2 0;...
           2 5 1 2 0;6 9 3 2 0; 6 13 1 2 0;6 14 2 2 0;6 13 1 2 0;5 11 1 2 1;...
           4 8 2 3 0;5 9 1 2 0; 5 10 2 3 0;6 14 4 2 0; 3 7 1 3 0; 4 9 1 3 0;...
           5 11 1 2 0; 11 12 2 2 0;9 11 1 3 0];  
molCompP = zeros(1,5);       
for ind=1:length(aminoacids)
    rep = length(find(seq==aminoacids(ind)));
    molCompP = molCompP + rep*amiComp(ind,:);
end
    
    