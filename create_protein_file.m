
AA = fastaread('protein\protein_all.fasta');

virusProtFinal = {};
for prot = 1%:length(AA)
    virusProt = {};
    header = AA(prot).Header;
    sequence = AA(prot).Sequence;
    idx = findstr(header,'|'); % find indices of header delimiter
    virusProt = [virusProt , header(1:(idx(1)-1))];
    for ii = 2:(length(idx)-1)
        virusProt = [virusProt , header(idx(ii)+1:idx(ii+1)-1)];
    end
    virusProt = [virusProt, header((idx(end)+1):end), AA(prot).Sequence];
    virusProtFinal = [virusProtFinal;virusProt]
end
    
        
        