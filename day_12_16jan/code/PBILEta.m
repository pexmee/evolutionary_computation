function PBILEta
nMonte = 20;
Problem = @Ackley;
Display = false;
GenLimit = 100;
Eta = [0.02, 0.1, 0.2];
MinCostEta1 = zeros(GenLimit+1, nMonte);
MinCostEta2 = zeros(GenLimit+1, nMonte);
MinCostEta3 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostEta1(:,i) = PBIL(Problem, Display, GenLimit, Eta(1));
    MinCostEta2(:,i) = PBIL(Problem, Display, GenLimit, Eta(2));
    MinCostEta3(:,i) = PBIL(Problem, Display, GenLimit, Eta(3));
end
MinCostEta1 = mean(MinCostEta1, 2);
MinCostEta2 = mean(MinCostEta2, 2);
MinCostEta3 = mean(MinCostEta3, 2);
close all
figure
plot(0:GenLimit,MinCostEta1,'r--', 0:GenLimit,MinCostEta2,'k:', 0:GenLimit,MinCostEta3,'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend(['Learning Rate = ',num2str(Eta(1))], ['Learning Rate = ',num2str(Eta(2))], ...
    ['Learning Rate = ',num2str(Eta(3))])