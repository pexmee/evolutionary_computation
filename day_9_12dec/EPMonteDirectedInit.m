function EPMonteDirectedInit(ProblemFunction, nMonte)
% Monte Carlo simulation to explore the effect of direction initialization using EP
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20;
end
DisplayFlag = false;
AdaptFlag = false;
PopSize = 10;
MaxGen = 50;
MinCostLessInit = zeros(MaxGen+1, nMonte);
MinCostMoreInit = zeros(MaxGen+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = round(sum(clock));
    MinCostLessInit(:, i) = EP(ProblemFunction, DisplayFlag, AdaptFlag, PopSize, PopSize, RandSeed, MaxGen);
    MinCostMoreInit(:, i) = EP(ProblemFunction, DisplayFlag, AdaptFlag, PopSize, 2*PopSize, RandSeed, MaxGen);
end
MinCostLessInit = mean(MinCostLessInit, 2);
MinCostMoreInit = mean(MinCostMoreInit, 2);
figure
plot((0 : MaxGen), MinCostLessInit, 'r--', (0 : MaxGen), MinCostMoreInit, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Without Extra Initial Individuals', 'With Extra Initial Individuals')