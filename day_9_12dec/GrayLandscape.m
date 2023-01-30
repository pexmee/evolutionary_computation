function GrayLandscape
% Create Figure 8.3
SetPlotOptions
close all
x = -4 : 5/15 : 1+eps;
f = inline('x.^4+5*x.^3+4*x.^2-4*x+1');
fArr = f(x);
fbinary = [fArr(1) fArr(2) fArr(4) fArr(3) fArr(8) fArr(7) fArr(5) fArr(6) ...
    fArr(16) fArr(15) fArr(13) fArr(14) fArr(9) fArr(10) fArr(12) fArr(11)];
figure; hold on; box on
subplot(2,1,1)
plot(x, fbinary, 'o') 
xlabel('binary-coded x')
ylabel('f(x)')
set(gca, 'XTick', x)
%set (gca, 'xticklabel', ['-4';'  ';'  ';'-3';'  ';'  ';'-2';'  ';'  ';'-1';'  ';'  ';' 0';'  ';'  ';' 1']);
set (gca, 'xticklabel', []);
set(gca, 'yticklabel', [])
set(gca, 'YTick', [])
%set(gca, 'xticklabel', ['0000'; '0001'; '0010'; '0011'; ...
%    '0100'; '0101'; '0110'; '0111'; ...
%    '1000'; '1001'; '1010'; '1011'; ...
%    '1100'; '1101'; '1110'; '1111'])
subplot(2,1,2)
plot(x, fArr, 'o')
xlabel('gray-coded x')
ylabel('f(x)')
%set (gca, 'xticklabel', ['-4';'  ';'  ';'-3';'  ';'  ';'-2';'  ';'  ';'-1';'  ';'  ';' 0';'  ';'  ';' 1']);
set (gca, 'xticklabel', []);
set(gca, 'XTick', x)
set(gca, 'yticklabel', [])
set(gca, 'YTick', [])