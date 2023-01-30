function MonteESmulambdaAdapt(ProblemFunction)
% Monte Carlo simulation to compare an ES with mutation rate adaptation and an ES without it
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
nMonte = 100;
Mu = 10;
Lambda = 20;
GenLimit = 100;
DisplayFlag = false;
PlusFlag = true;
MinCostNoAdapt = zeros(GenLimit+1, nMonte);
MinCostYesAdapt = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostNoAdapt(:,i) = ES(ProblemFunction, DisplayFlag, PlusFlag, 0, Mu, Lambda, GenLimit); % no adapt
    MinCostYesAdapt(:,i) = ES(ProblemFunction, DisplayFlag, PlusFlag, 2, Mu, Lambda, GenLimit); % yes adapt
end
MinCostNoAdapt = mean(MinCostNoAdapt, 2);
MinCostYesAdapt = mean(MinCostYesAdapt, 2);
figure
plot(0:GenLimit, MinCostNoAdapt, 'r--', 0:GenLimit, MinCostYesAdapt, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Without Mutation Rate Adaptation', 'With Mutation Rate Adaptation')