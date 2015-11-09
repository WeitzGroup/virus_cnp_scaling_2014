function molCompDNA = get_DNA_mol_form(phageName)
        
    % creates the name phage_DNA.txt and gets the sequence in that file
    filename = [phageName '_DNA.txt'];
    
    nucleobase = fastaread(filename);
    molCompDNA = DNA_mol_form(nucleobase.Sequence);    
    DNAStruct =2*length(nucleobase.Sequence)*[5 0 1];
    molCompDNA = molCompDNA + DNAStruct;  
end



