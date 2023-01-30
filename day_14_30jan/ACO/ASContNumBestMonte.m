function ASContNumBestMonte(nMonte)
% Monte Carlo simulation of ant system code to explore the effect of the number of pheromone contributors
% INPUT: nMonte = number of Monte Carlo simulations
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20;
end
NumIntervals = 20;
PopSize = 40;
GenLimit = 100;
tauMin = 0;
tauMax = inf;
DisplayFlag = false;
MinCostNumBest40 = zeros(GenLimit+1, nMonte);
MinCostNumBest4 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostNumBest40(:,i) = ASCont(@Ackley, DisplayFlag, NumIntervals, PopSize, GenLimit, tauMin, tauMax, 40);
    MinCostNumBest4(:,i) = ASCont(@Ackley, DisplayFlag, NumIntervals, PopSize, GenLimit, tauMin, tauMax, 4);
end
MinCostNumBest40 = mean(MinCostNumBest40, 2);
MinCostNumBest4 = mean(MinCostNumBest4, 2);
SetPlotOptions
figure
plot(0:GenLimit, MinCostNumBest40, 'r--', 0:GenLimit, MinCostNumBest4, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('40 Pheromone Contributors', '4 Pheromone Contributors')