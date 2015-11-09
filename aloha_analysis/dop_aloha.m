%botid, date, time, press, dop, 
%#, mmddyy, hhmmss, dbar, umol/kg, 

% Read in data
x=tdfread('dop_data_justnums.txt');   
% Clean up dates
for i=1:length(x.Date),
  if (x.Date(i)>0)
    x.year(i)=mod(mod(x.Date(i),10000),100);
    tmpd = (x.Date(i)-x.year(i))/100;  % remove last two digits
    x.day(i)=rem(tmpd,100);
    tmpd = (tmpd-x.day(i))/100;
    x.month(i)=tmpd;
    if (x.year(i)>50)
      x.year(i)=1900+x.year(i);
    else 
      x.year(i)=2000+x.year(i);
    end
    x.serial_date(i)=datenum(x.year(i),x.month(i),x.day(i));
    x.date_valid(i)=1;
  else
    x.date_valid(i)=0;
    x.month(i)=-1;
    x.day(i)=-1;
    x.year(i)=-1;
  end
end
tmpi=find(x.date_valid);
plot(x.serial_date(tmpi),x.DOP(tmpi),'ko');




