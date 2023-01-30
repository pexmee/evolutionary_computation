function MonteCGABinaryAlpha
nMonte = 50;
GenLimit = 50;
DisplayFlag = false;
MinCostAlpha1 = zeros(GenLimit+1, nMonte);
MinCostAlpha2 = zeros(GenLimit+1, nMonte);
MinCostAlpha3 = zeros(GenLimit+1, nMonte);
MinCostAlpha4 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    MinCostAlpha1(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 0.1); 
    MinCostAlpha2(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 0.01); 
    MinCostAlpha3(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 0.001); 
    MinCostAlpha4(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 0.0001); 
end
MinCostAlpha1 = mean(MinCostAlpha1, 2);
MinCostAlpha2 = mean(MinCostAlpha2, 2);
MinCostAlpha3 = mean(MinCostAlpha3, 2);
MinCostAlpha4 = mean(MinCostAlpha4, 2);
figure, hold on
plot(0:GenLimit, MinCostAlpha1, 'k-.')
plot(0:GenLimit, MinCostAlpha2, 'b-')
plot(0:GenLimit, MinCostAlpha3, 'r:')
plot(0:GenLimit, MinCostAlpha4, 'g-.')
xlabel('Generation')
ylabel('Minimum Cost')
legend('alpha = 0.1', 'alpha = 0.01', 'alpha = 0.001', 'alpha = 0.0001')
box on