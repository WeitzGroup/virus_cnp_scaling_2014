function molForm = DNA_mol_form(seq)
%molComP = [C N P]
%molecular composition of DNA whose bases are specified in seq
% inlucing the complementary bases
base = 'ATCG';
baseComp = [5 5 0; 5 2 0; 4 3 0 ;5 5 0 ];
basePairComp = [ 5 2 0; 5 5 0; 5 5 0; 4 3 0]; %complementary base composition 'TAGC'
molForm = zeros(1,3);       
for ind=1:length(base)
    rep = length(find(seq==base(ind)));
    molForm = molForm + rep*(baseComp(ind,:) + basePairComp(ind,:));
end

