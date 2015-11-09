function molCompCap = get_cap_mol_comp(rdata,phageName)
% in: raw excel data, name of phage
% out: 
    capComp = get_capComp(phageName);
    [typesProt,~] =size(capComp);
    molCompCap = zeros(1,5);
    matlabDir = dir('*.txt');
    
    for j = 1:typesProt
        %capComp{j,1}
        sequence =  get_sequence(phageName,capComp{j,1});
       
        nProt = capComp{j,2};
        molCompProt = protein_mol_form(sequence);  %[C H N O S]
        molCompCap = molCompCap + nProt*molCompProt;
    end
    
    
    function cComp = get_capComp(phageName)
    % in: str phageName     
    % out: cell 2by2 first col are the names of the protein, second the number of each one.  
    % uses excel file with specific location of names and types of protein.
        [row,column]= find(strcmp(rdata,phageName));
        numberProt = rdata{row, column+1};
        cComp =rdata(row:(row+numberProt-1),3:4);
    end    
    
    
    function sequence =  get_sequence(phageName,protein)
    %  creates phage_protein.txt and gets the sequence of AA in that file
    % in: str; phageName. str; protein
    % out: creates the name 
        filename = [phageName '_' protein '.txt'];
        for i = 1: length(matlabDir);
            if strcmp(matlabDir(i).name, filename);
                AA = fastaread(matlabDir(i).name);
                %%AA.Header
                sequence = AA.Sequence;
                break
            end
        end
    end 

end