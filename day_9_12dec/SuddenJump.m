function SuddenJump
gen = 0 : 100;
cost = zeros(size(gen));
cost(1) = 10;
for i = 2 : 70
    sigma = exp(-i/10);
    cost(i) = cost(i-1) - abs(sigma * randn);
end
jump = i + 1;
cost(jump) = cost(jump-1) - 1;
for i = jump+1 : length(gen)
    sigma = exp(-i/15);
    cost(i) = cost(i-1) - abs(sigma * randn);
end
SetPlotOptions
plot(gen, cost)
xlabel('Generation')
ylabel('Cost')
