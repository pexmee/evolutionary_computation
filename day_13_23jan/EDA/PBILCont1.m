function PBILCont1

sigmax = 5;
meanx = 3;
x = -3*sigmax : 0.1 : 3*sigmax;
pdfx1 = 1 / sqrt(2 * pi) / sigmax * exp(-(x - meanx).^2 / 2 / sigmax^2);

sigmax2 = 3;
meanx2 = -3;
pdfx2 = 1 / sqrt(2 * pi) / sigmax2 * exp(-(x - meanx2).^2 / 2 / sigmax2^2);

SetPlotOptions
plot(x,pdfx1,'r--', x,pdfx2,'b-');
minx = -2 * sigmax;
maxx = 2 * sigmax;
set(gca, 'XTick', [minx, maxx])
set(gca, 'xticklabel', ['min(x)'; 'max(x)'])
set(gca, 'YTick', [])
grid
text(4, 0.085, 'Initial pdf')
text(-1, 0.125, 'Final pdf')