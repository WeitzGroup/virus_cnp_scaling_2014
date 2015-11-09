protein_char_all    % script to calculate dC, dN, and respective errors
radius_vs_nbp_data  %script to calculate the fit for fillfit
fillC =roundn(fillfit,-2); % fraction of the capsid fill with DNA rounded to two decimals
r = 10:50;        %capsid radius
nbp = 880:180000;    %number of base pairs
Vbp = (1)^2*(0.34)*pi;   %volume of a base pair

mCbp = [19.5 7.5 2];   % C:N:P molecular formula of 2 nucleotides (A + T + C + G)/2 + [10 0 2]
NP_bp = mCbp(2)/mCbp(3); 
CN_bp = mCbp(1)/mCbp(2);
CP_bp = mCbp(1)/mCbp(3);
%CN_Pr = protMolComp(1)/protMolComp(2);

%% nbp as a function of outer radius in nm
nbpR = @(x)  4/3*pi*(x-h).^3./Vbp.*fillC; 
nbp_r =nbpR(r) ; 

%% theory, molecular content DNA as a function of n_bp and r
thDNA.C.r = mCbp(1)*nbp_r;
thDNA.N.r = mCbp(2)*nbp_r;
thDNA.P.r = mCbp(3)*nbp_r;

thDNA.C.bp = mCbp(1)*nbp;
thDNA.N.bp = mCbp(2)*nbp;
thDNA.P.bp = mCbp(3)*nbp;

errh = 1;
errCbp = 0.5;
errNbp = 1.5;

errDNA.C.r = sqrt( mCbp(1)^2*(r - h).^6*errC^2 + nbp_r.^2*errCbp^2 ...
                + (mCbp(1)*4/3*pi/Vbp*fillC)^2*9*(r-h).^4*errh^2); 
errDNA.N.r = sqrt( mCbp(2)^2*(r - h).^6*errC^2 +  ...
                + (mCbp(2)*4/3*pi/Vbp*fillC)^2*9*(r-h).^4*errh^2); 
errDNA.P.r = sqrt( mCbp(3)^2*(r - h).^6*errC^2  ...
                + (mCbp(3)*4/3*pi/Vbp*fillC)^2*9*(r-h).^4*errh^2); 

errDNA.C.bp = errCbp*nbp;
errDNA.N.bp = errNbp*nbp;


%% theory capsid

r_bp =@(x) ((x*Vbp*3)/(4*pi*fillC)).^(1/3)+ h;

capsid_carbon = @(x) 4*pi*dC*(x.^3 - (x-h).^3)/3;
capsid_nitrogen = @(x) 4*pi*dN*(x.^3 - (x-h).^3)/3;

thCapsid.C.r = capsid_carbon(r);
thCapsid.N.r = capsid_nitrogen(r);

thCapsid.C.bp =capsid_carbon(r_bp(nbp)); 
thCapsid.N.bp = capsid_nitrogen(r_bp(nbp));

%Errors
errCapsid.C.r = sqrt( (thCapsid.C.r/dC*errdC).^2 + ...
                (4*pi*dC*(3*r.^2 - 6*h*r.^2 + 3*h^2)/3*errh).^2);
errCapsid.N.r = sqrt( (thCapsid.C.r/dN*errdN).^2 + ...
                (4*pi*dN*(3*r.^2 - 6*h*r.^2 + 3*h^2)/3*errh).^2); 
 
drdV_F= 1/3*(nbp*3/4/pi).^(1/3)*(Vbp/fillC)^(-2/3);
c = fillC*4*pi/3/Vbp;            
errCapsid.C.bp =  sqrt( (thCapsid.C.bp/dC*errdC).^2 + ...
                (dC*4*pi*r_bp(nbp).^2*errh).^2 + (dC*(4*pi*r_bp(nbp).^2.*drdV_F - nbp)*4*pi/c^2*errC).^2) ;
errCapsid.N.bp =  sqrt( (thCapsid.N.bp/dN*errdN).^2 + ...
                (dN*4*pi*r_bp(nbp).^2*errh).^2 + (dN*(4*pi*r_bp(nbp).^2.*drdV_F - nbp)*4*pi/c^2*errC).^2) ;
            
%% molecular composition of phage
thPhage.C.r = thCapsid.C.r +thDNA.C.r;
thPhage.N.r = thCapsid.N.r +thDNA.N.r;
thPhage.P.r = thDNA.P.r;

errPhage.C.r = sqrt(errCapsid.C.r.^2 +errDNA.C.r.^2);
errPhage.N.r = sqrt(errCapsid.N.r.^2 +errDNA.N.r.^2);
errPhage.P.r = errDNA.P.r;

thPhage.C.bp = thCapsid.C.bp + thDNA.C.bp;
thPhage.N.bp = thCapsid.N.bp + thDNA.N.bp ;
thPhage.P.bp=  thDNA.P.bp;

errPhage.C.bp = sqrt(errCapsid.C.bp.^2 +errDNA.C.bp.^2);
errPhage.N.bp = sqrt(errCapsid.N.bp.^2 +errDNA.N.bp.^2);
