function amiFreq = amino_freq(seq)
%molComP = [C H N O S]
aminoacids = 'ACDEFGHIKLMNPQRSTVWY';  
amiFreq = zeros(1,length(aminoacids));
for ind=1:length(aminoacids)
    amiFreq(ind) = length(find(seq==aminoacids(ind)));
end
end