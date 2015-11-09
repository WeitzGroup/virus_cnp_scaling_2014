% Read in digitized data
dip=tdfread('dip_measurements.csv',',');
vir=tdfread('vcounts_matteson.csv',',');
% Correct for the missing data on day 2 for DIP
tmpday=dip.Day;
tmpday(3:16)=tmpday(2:15);
tmpday(2)=vir.Day(2);
tmpdip=dip.DRP;
tmpdip(3:16)=tmpdip(2:15);
tmpdip(2)=-1;
dip.Day=vir.Day;
dip.DRP=tmpdip*1000;
% Reference
vir.dop_conc_high=vir.Virus/(5*10^9);    % If it is 5*10^9 it is 1nm
vir.dop_conc_low=vir.Virus/(5*10^9)*0.3; % If it is 5*10^9 it is 0.3 nm
vir.dop_highratio = vir.dop_conc_high./dip.DRP;
vir.dop_lowratio = vir.dop_conc_low./dip.DRP;

