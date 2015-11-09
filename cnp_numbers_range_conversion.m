
avo = 6.02214179e23;
atmW = [12.01  14.01 30.97]; % [C N P] atomic weigth g/mol
viralDen = 10^11; %viruses per L
radius_vs_nbp_data
fillC = roundn(fillfit,-2);
h = 2.5;
Vbp = 1*0.34*pi;
r_bp =@(x) ((x*Vbp*3)/(4*pi*fillC)).^(1/3)+ h;

%% Virus-induced transformation of elemental ratios(stoichiometric mismatch)
med4.c=[45.8 60.9];
med4.n=[9.4 9.6];
med4.p=[0.98 0.34]; %[replete limited]
v.bp=50000;
v.r=35; % nm
v.burst=40;

[v.cmol,v.nmol,v.pmol] = cnp_head_r(v.r);
v.c =  v.cmol.*atmW(1)./avo.*10^15
v.n =  v.nmol.*atmW(2)./avo.*10^15
v.p = v.pmol*atmW(3)./avo.*10^15
lys.c = med4.c(1) -v.burst*v.c;
lys.n = med4.n(1) -v.burst*v.n;
lys.p = med4.p(1) -v.burst*v.p;

med4Ratiofg = [med4.c(1) med4.n(1) med4.p(1)]./med4.p(1)
vRatiofg = [v.c v.n v.p]./v.p
lysRatiofg = [lys.c lys.n lys.p]./lys.p

med4RatioMol = [med4.c(1)/atmW(1) med4.n(1)/atmW(2) med4.p(1)/atmW(3)]/(med4.p(1)/atmW(3))
vRatioMol = [v.cmol v.nmol v.pmol]./v.pmol
lysRatioMol = [lys.c/atmW(1) lys.n/atmW(2) lys.p/atmW(3)]/(lys.p/atmW(3))

fracVirPLimited = v.burst*v.p/med4.p(2)
[~,~,v.pmolOld] = cnp_head_r_old(v.r);
fracVirPLimitedOld = v.burst*v.pmolOld*atmW(3)./avo.*10^15/med4.p(2)
%% figure3 bars

he = 36.259; % height of full bar in px;

vBluebar = he*v.c*v.burst/med4.c(1)
vGreenbar = he*v.n*v.burst/med4.n(1)
vRedbar = he*v.p*v.burst/med4.p(1)

lysBluebar = he*lys.c/med4.c(1)
lysGreebar = he*lys.n/med4.n(1)
lysRedbar  = he*lys.p/med4.p(1)

%% cnp per particle in fg for range determined by dna length in data (776)
ri =r_bp(5475)
rf = r_bp(358663)

[vc,vn,vp] = cnp_head_r([ri rf]);

cfg =  vc.*atmW(1)./avo.*10^15
nfg =  vn.*atmW(2)./avo.*10^15
pfg = vp.*atmW(3)./avo.*10^15

%% ranges of capsid diameters
load DNA_sequences
for ind= 1 :length(nucl);
    dsLen(ind)  = length(nucl(ind).Sequence);
end

dsLen = sort(dsLen);
r_inf = r_bp(dsLen);
medianDiam = quantile(r_inf,.5)*2
range95diam = [quantile(r_inf,.025) quantile(r_inf,.975)]*2

%% c n p in fg for 50-70 nm range
[vc,vn,vp] = cnp_head_r([25 35]);
cfg2 =  vc.*atmW(1)./avo.*10^15
nfg2 =  vn.*atmW(2)./avo.*10^15
pfg2 = vp.*atmW(3)./avo.*10^15

%% P contribution in nM for densities 10^9-10^10

viralDen = [ 10^9 10^10];
[vc,vn,vp] = cnp_head_r([25 35]);

PnmolL1 = [vp(1)/avo*viralDen(1) vp(2)/avo*viralDen(2)]*10^9
CnomlL  = [ vc(1)*viralDen(1) vc(2)*viralDen(2)]*10^9/avo
NnomlL  = [ vn(1)*viralDen(1) vn(2)*viralDen(2)]*10^9/avo
%% P contribution in nM for densities 10^9-10^10
viralDen = [ 10^8 10^11];
[vc,vn,vp] = cnp_head_r([25 35]);

PnmolL1 = [vp(1)/avo*viralDen(1) vp(2)/avo*viralDen(2)]*10^9
CnomlL  = [ vc(1)*viralDen(1) vc(2)*viralDen(2)]*10^9/avo
NnomlL  = [ vn(1)*viralDen(1) vn(2)*viralDen(2)]*10^9/avo
%% Sargasso contribution of viruses to DOP
lowLimit = 0.079/100
highLimit = 2.4/50

%% HOT
viralDen3 = [5]*10^9;
[vc,vn,vp] = cnp_head_r([25 35]);
PnmolL2 = [vp(1)/avo*viralDen3(1) vp(2)/avo*viralDen3(1)]*10^9

contribVir = [0.4/230 1.2/150] % percentage contribution to total DOP pool
%% Brum ALOHA
viralDen4 = [8 11.3]*10^9;
[vc,vn,vp] = cnp_head_r([25 35]);
PnmolL3 = [vp(1)/avo*viralDen4(2) vp(2)/avo*viralDen4(1)]*10^9
contribVir = [0.9 1.9]./224 % percentage contribution to total DOP pool
%% southern pacific ocean
viralDen2 = [1 1.2]*10^11;
[vc,vn,vp] = cnp_head_r([25 35]);
PnmolL2 = [vp(1)/avo*viralDen2(1) vp(2)/avo*viralDen2(2)]*10^9

%% viral density to get 5nM DOP
[vc,vn,vp] = cnp_head_r([35]);

vD = 5*10^-9/vp*avo


