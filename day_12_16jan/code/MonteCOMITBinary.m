function MonteCOMITBinary(ProblemFunction, EntropyFlag)
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @AckleyDisc;
end
if ~exist('EntropyFlag', 'var') || isempty(EntropyFlag)
    EntropyFlag = true;
end
nMonte = 20;
GenLimit = 20;
DisplayFlag = false;
MinCostA = zeros(GenLimit+1, nMonte);
MinCostB = zeros(GenLimit+1, nMonte);
MinCostC = zeros(GenLimit+1, nMonte);
BestCostA = zeros(GenLimit+1, 1);
BestCostB = zeros(GenLimit+1, 1);
BestCostC = zeros(GenLimit+1, 1);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    MinCostA(:, i) = MIMICBinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, 40, true, false, EntropyFlag); 
    MinCostB(:, i) = MIMICBinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, 40, false, true, EntropyFlag); 
    MinCostC(:, i) = MIMICBinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, 40, false, false, EntropyFlag); 
end
for i = 1 : GenLimit+1
    BestCostA(i) = min(MinCostA(i,:));
    BestCostB(i) = min(MinCostB(i,:));
    BestCostC(i) = min(MinCostC(i,:));
end
MinCostA = mean(MinCostA, 2);
MinCostB = mean(MinCostB, 2);
MinCostC = mean(MinCostC, 2);
% Plot average of best cost at each generation
figure, hold on
plot(0:GenLimit, MinCostA, 'k-.')
plot(0:GenLimit, MinCostB, 'r:')
plot(0:GenLimit, MinCostC, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
title('Average Monte Carlo Results')
if EntropyFlag
    legend('First Order', 'Random Permutation', 'MIMIC')
else
    legend('First Order', 'Random Permutation', 'COMIT')
end
% Plot best of best cost at each generation
figure, hold on
plot(0:GenLimit, BestCostA, 'k-.')
plot(0:GenLimit, BestCostB, 'r:')
plot(0:GenLimit, BestCostC, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
title('Best Monte Carlo Results')
if EntropyFlag
    legend('First Order', 'Random Permutation', 'MIMIC')
else
    legend('First Order', 'Random Permutation', 'COMIT')
end