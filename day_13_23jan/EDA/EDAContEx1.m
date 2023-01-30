function EDAContEx1

x = 0 : 0.1 : 1;
y = 1.5 - x;
x = [-0.5 0 x 1 1.5];
y = [ 0 0 y 0 0];
SetPlotOptions
close all
figure
plot(x, y)
xlabel('x(i)')
ylabel('PDF of x(i)')
axis([-0.5 1.5 0 2])