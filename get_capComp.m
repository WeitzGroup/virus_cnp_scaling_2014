function cComp = get_capComp(rdata,phageName)
    % in: str phageName     
    % out: cell 2by2 first col are the names of the protein, second the number of each one.  
    % uses excel file with specific location of names and types of protein.
        [row,column]= find(strcmp(rdata,phageName));
        numberProt = rdata{row, column+1};
        cComp =rdata(row:(row+numberProt-1),3:4);
    end    