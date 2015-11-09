Id 
Sample Id
decy 
Decimal Year
Deployed
decd 
Decimal Day
Deployed
qc 
Niskin
Quality Flag
Depth 
Depth (m)
PO4 
Phosphate 
(µmol/kg)
TDP 
Total Dissolved 
Phosphorus 
(nmol/kg)
 	Data Quality flags:	-9.99 = No data
 	0 = Less than detection limit
 	Bottle Quality flags:	1 = Normal Bottle Fire
 	2 = Bottle Misfire
 	3 = Suspect Bottle
 	All positions are in decimal degrees.
 	All times are GMT.

Matlab file has 7 columns
ID
Decimal year
Decimal day
Quality
Depth
PO4 (um)
TDP (nmol)

% Depth against DOP
% suggests range from 10-100
plot(data(tmpi,5),data(tmpi,7)/1000-data(tmpi,6),'ko');

