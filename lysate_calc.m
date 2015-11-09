function lysate = lysate_calc(iso,v,x)
% function lysate = lysate_calc(iso,v,x)
% 
% iso is the isolate (structure c, n and p, with ranges)
% v is the virus, with radius, bp, and burst size)
% x is the range of conditions
%
% lysate is the stuff released, with c and n and p in ranges
for i=1:length(x),
  r=x(i);
  dc=(iso.c(2)-iso.c(1));
  dn=(iso.n(2)-iso.n(1));
  dp=(iso.p(2)-iso.p(1));
  c=iso.c(1)+r*dc;
  n=iso.n(1)+r*dn;
  p=iso.p(1)+r*dp;
  lysate.c(i)=c-v.burst*v.c;
  lysate.n(i)=n-v.burst*v.n;
  lysate.p(i)=p-v.burst*v.p;
  lysate.isoc(i)=c;
  lysate.ison(i)=n;
  lysate.isop(i)=p;
end
lysate.level=x;
end
