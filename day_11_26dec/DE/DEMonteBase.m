function DEMonteBase(nMonte)
% Compare DE performance using a random individual, the current individual, or the best individual as the base vector
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20; % number of Monte Carlo simulations
end
Problem = @Ackley;
Display = false;
GenLimit = 80;
LFlag = false; % Use the "/bin" version of DE
NumDiff = 1; % Use 1 difference vector
MinCostBase1 = zeros(GenLimit+1, nMonte);
MinCostBase2 = zeros(GenLimit+1, nMonte);
MinCostBase3 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostBase1(:,i) = DE(Problem, Display, GenLimit, LFlag, 1, NumDiff);
    MinCostBase2(:,i) = DE(Problem, Display, GenLimit, LFlag, 2, NumDiff);
    MinCostBase3(:,i) = DE(Problem, Display, GenLimit, LFlag, 3, NumDiff);
end
MinCostBase1 = mean(MinCostBase1, 2);
MinCostBase2 = mean(MinCostBase2, 2);
MinCostBase3 = mean(MinCostBase3, 2);
SetPlotOptions
figure
plot(0:GenLimit,MinCostBase2,'r--', 0:GenLimit,MinCostBase3,'k:', 0:GenLimit,MinCostBase1,'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Random base vector', 'Current base vector', 'Best base vector')
