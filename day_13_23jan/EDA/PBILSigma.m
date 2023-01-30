function PBILSigma
nMonte = 50;
Problem = @Ackley;
Display = false;
GenLimit = 200;
Eta = 0.1;
UpdateCount = 5;
Sigma = [1/50, 1/50; 1/10, 1/10; 1/10, 1/50];
MinCostSigma1 = zeros(GenLimit+1, nMonte);
MinCostSigma2 = zeros(GenLimit+1, nMonte);
MinCostSigma3 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostSigma1(:,i) = PBIL(Problem, Display, GenLimit, Eta, UpdateCount, Sigma(1,:));
    MinCostSigma2(:,i) = PBIL(Problem, Display, GenLimit, Eta, UpdateCount, Sigma(2,:));
    MinCostSigma3(:,i) = PBIL(Problem, Display, GenLimit, Eta, UpdateCount, Sigma(3,:));
end
MinCostSigma1 = mean(MinCostSigma1, 2);
MinCostSigma2 = mean(MinCostSigma2, 2);
MinCostSigma3 = mean(MinCostSigma3, 2);
close all
figure
plot(0:GenLimit,MinCostSigma1,'r--', 0:GenLimit,MinCostSigma2,'k:', 0:GenLimit,MinCostSigma3,'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend(['k_0 = ', num2str(Sigma(1,1)), ', k_f = ', num2str(Sigma(1,2))], ...
    ['k_0 = ', num2str(Sigma(2,1)), ', k_f = ', num2str(Sigma(2,2))], ...
    ['k_0 = ', num2str(Sigma(3,1)), ', k_f = ', num2str(Sigma(3,2))]);