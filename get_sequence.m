function sequence =  get_sequence(matlabDir,phageName,protein)
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