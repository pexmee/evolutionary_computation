function MonteESmulambda(ProblemFunction)
% Monte Carlo simulation to compare a (mu+lambda)-ES with a (mu,lambda)-ES
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
nMonte = 100;
Mu = 10;
Lambda = 20;
GenLimit = 100;
DisplayFlag = false;
AdaptFlag = false;
MinCostPlus = zeros(GenLimit+1, nMonte);
MinCostComma = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostPlus(:,i) = ES(ProblemFunction, DisplayFlag, true, AdaptFlag, Mu, Lambda, GenLimit); % mu+lambda
    MinCostComma(:,i) = ES(ProblemFunction, DisplayFlag, false, AdaptFlag, Mu, Lambda, GenLimit); % mu,lambda
end
MinCostPlus = mean(MinCostPlus, 2);
MinCostComma = mean(MinCostComma, 2);
figure
plot(0:GenLimit, MinCostPlus, 'r--', 0:GenLimit, MinCostComma, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('(mu+lambda)-ES', '(mu,lambda)-ES')