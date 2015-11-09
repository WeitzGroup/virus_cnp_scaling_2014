function printpdf(filename)
% creates a pdf and a eps of the current figure

print([filename '.eps'],'-depsc2');
system(['epstopdf ',[filename '.eps']]);
delete([filename '.eps'])
