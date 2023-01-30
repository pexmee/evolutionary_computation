function PBILUpdateCount
nMonte = 50;
Problem = @Ackley;
Display = false;
GenLimit = 200;
Eta = 0.1;
UpdateCount = [1, 10, 20];
MinCostUpdate1 = zeros(GenLimit+1, nMonte);
MinCostUpdate2 = zeros(GenLimit+1, nMonte);
MinCostUpdate3 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostUpdate1(:,i) = PBIL(Problem, Display, GenLimit, Eta, UpdateCount(1));
    MinCostUpdate2(:,i) = PBIL(Problem, Display, GenLimit, Eta, UpdateCount(2));
    MinCostUpdate3(:,i) = PBIL(Problem, Display, GenLimit, Eta, UpdateCount(3));
end
MinCostUpdate1 = mean(MinCostUpdate1, 2);
MinCostUpdate2 = mean(MinCostUpdate2, 2);
MinCostUpdate3 = mean(MinCostUpdate3, 2);
close all
figure
plot(0:GenLimit,MinCostUpdate1,'r--', 0:GenLimit,MinCostUpdate2,'k:', 0:GenLimit,MinCostUpdate3,'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend(['# Best = # Worst = ',num2str(UpdateCount(1))], ['# Best = # Worst = ',num2str(UpdateCount(2))], ...
    ['# Best = # Worst = ',num2str(UpdateCount(3))]);