function MonteES1plus1(ProblemFunction, nMonte)
% Monte Carlo simulation to compare an ES with standard deviation adaptation with an ES without standard deviation adaptation
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 100;
end
Mu = 1;
Lambda = 1;
GenLimit = 500;
DisplayFlag = false;
PlusFlag = true;
MinCostNoAdapt = zeros(GenLimit+1, nMonte);
MinCostYesAdapt = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Monte Carlo simulation ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostNoAdapt(:,i) = ES(ProblemFunction, DisplayFlag, PlusFlag, false, Mu, Lambda, GenLimit); % no sigma adaptation
    MinCostYesAdapt(:,i) = ES(ProblemFunction, DisplayFlag, PlusFlag, true, Mu, Lambda, GenLimit); % yes sigma adaptation
end
MinCostNoAdapt = mean(MinCostNoAdapt, 2);
MinCostYesAdapt = mean(MinCostYesAdapt, 2);
figure
plot(0:GenLimit, MinCostNoAdapt, 'r--', 0:GenLimit, MinCostYesAdapt, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Without Standard Deviation Adaptation', 'With Standard Deviation Adaptation')