function MonteUMDABinary
nMonte = 50;
GenLimit = 50;
DisplayFlag = false;
MinCostApct = zeros(GenLimit+1, nMonte);
MinCostBpct = zeros(GenLimit+1, nMonte);
MinCostCpct = zeros(GenLimit+1, nMonte);
MinCostDpct = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    MinCostApct(:, i) = UMDABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 10); 
    MinCostBpct(:, i) = UMDABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 40); 
    MinCostCpct(:, i) = UMDABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 70); 
    MinCostDpct(:, i) = UMDABinary(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, 2); 
end
MinCostApct = mean(MinCostApct, 2);
MinCostBpct = mean(MinCostBpct, 2);
MinCostCpct = mean(MinCostCpct, 2);
MinCostDpct = mean(MinCostDpct, 2);
figure, hold on
plot(0:GenLimit, MinCostDpct, 'm-.')
plot(0:GenLimit, MinCostCpct, 'b-')
plot(0:GenLimit, MinCostBpct, 'r:')
plot(0:GenLimit, MinCostApct, 'k-.')
xlabel('Generation')
ylabel('Minimum Cost')
legend('M = 2% of Population', 'M = 70% of Population', 'M = 40% of Population', 'M = 10% of Population')