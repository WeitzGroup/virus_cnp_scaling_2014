phageCol.name = {'T4','N4','Syn5'};
phageCol.r = [46.3 34.8 30];

capsidCol.C = [4226432 1137264 702398 ];
capsidCol.N = [1120272 304304  196376 ];

phageCol.C = [8713556 2998609 1715042 ];
phageCol.N = [2678916 956895 576459 ];
phageCol.P = [337800 140306 92428 ];


DNAcol.C = [3318381 1374084 898876];								 
DNAcol.N = [1241919 520047 348902];
DNAcol.P = [337800 140306 92428];

phageCol.length = DNAcol.P/2;

headCol.C = capsidCol.C + DNAcol.C;
headCol.N = capsidCol.N + DNAcol.N;
headCol.P = DNAcol.P;



