% A script to estimate the fraction of DOP bound
% in virus populations, using
%
% 1. Virus population ranges
% 2. DOP measurements

% Calibration
vhead.mweightP = 30.97; % g/mol
vhead.fgP =  [0.0025 0.0074]; % fg
vhead.nMP = vhead.fgP/vhead.mweightP*10^-15*10^9; % P in units of nM

% Data background
vfrac{1}.v = [6 12]*10^9;
vfrac{1}.dop = [36 80];
vfrac{1}.site = 'BATS - late summer';
vfrac{2}.v = [1 3]*10^9;
vfrac{2}.dop = [36 80];
vfrac{2}.site = 'BATS - stratified';
vfrac{3}.v = [4.5 5.5]*10^9;
vfrac{3}.dop = [210 250];
vfrac{3}.site = 'HOT - Culley';
vfrac{4}.v = [8 11.3]*10^9;
vfrac{4}.dop = [ 224-46 224+46];
vfrac{4}.site = 'HOT - Brum';
vfrac{5}.v = [17 120]*10^9;
vfrac{5}.dop = [150 225];
vfrac{5}.site = 'S. Pacific - Wilhelm';
vfrac{6}.v = [6.1 26]*10^9;
vfrac{6}.dop = [150 225];
vfrac{6}.site = 'S. Pacific - Brussaard, 2009';
vfrac{7}.v = [0.1 7.6]*10^9;
vfrac{7}.dop = [150 225];
vfrac{7}.site = 'S. Pacific - Brussaard, 2012';


% Fractional calculations
for i=1:length(vfrac),
  vfrac{i}.dop_v=[vfrac{i}.v(1)*vhead.nMP(1) vfrac{i}.v(2)*vhead.nMP(2)];
  vfrac{i}.frac_dop_v = [vfrac{i}.dop_v(1)/vfrac{i}.dop(2) vfrac{i}.dop_v(2)/vfrac{i}.dop(1)];
end
