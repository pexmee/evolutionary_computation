function MonteEAVarGA(Problem, nMonte)
% Explore the effect of gray-coding vs. binary-coding in a GA
if ~exist('Problem', 'var') || isempty(Problem)
    Problem = @WorstCaseProblem;
end
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 50;
end
GenLimit = 50;
DisplayFlag = false;
AveCostYesGray = zeros(GenLimit+1, nMonte);
AveCostNoGray = zeros(GenLimit+1, nMonte);
MinCostYesGray = zeros(GenLimit+1, nMonte);
MinCostNoGray = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    [MinCostYesGray(:, i), AveCostYesGray(:, i)] = GA(Problem, DisplayFlag, RandSeed, GenLimit, true, 0); 
    [MinCostNoGray(:, i), AveCostNoGray(:, i)] = GA(Problem, DisplayFlag, RandSeed, GenLimit, false, 0); 
end
AveCostYesGray = mean(AveCostYesGray, 2);
AveCostNoGray = mean(AveCostNoGray, 2);
MinCostYesGray = mean(MinCostYesGray, 2);
MinCostNoGray = mean(MinCostNoGray, 2);
figure, hold on, box on
plot(0:GenLimit, AveCostYesGray, 'b--')
plot(0:GenLimit, AveCostNoGray, 'r-')
xlabel('Generation')
ylabel('Average Cost')
legend('Gray Coding', 'Binary Coding')
figure, hold on, box on
plot(0:GenLimit, MinCostYesGray, 'b--')
plot(0:GenLimit, MinCostNoGray, 'r-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Gray Coding', 'Binary Coding')