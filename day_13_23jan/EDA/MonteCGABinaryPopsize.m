function MonteCGABinaryPopsize
nMonte = 50;
GenLimit = 50;
alpha = 0.001;
DisplayFlag = false;
MinCost1 = zeros(GenLimit+1, nMonte);
MinCost2 = zeros(GenLimit+1, nMonte);
MinCost3 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    MinCost1(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, alpha, 2); 
    MinCost2(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, alpha, 5); 
    MinCost3(:, i) = cGABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, alpha, 20); 
end
MinCost1 = mean(MinCost1, 2);
MinCost2 = mean(MinCost2, 2);
MinCost3 = mean(MinCost3, 2);
figure, hold on
plot(0:GenLimit, MinCost1, 'k-.')
plot(0:GenLimit, MinCost2, 'b-')
plot(0:GenLimit, MinCost3, 'r:')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Pop. Size = 2', 'Pop. Size = 5', 'Pop. Size = 20')