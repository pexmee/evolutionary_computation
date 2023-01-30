function ASContMonte(nMonte)
% Monte Carlo simulation of ant system code to explore the effect of the number of pheromone bins
% INPUT: nMonte = number of Monte Carlo simulations
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20;
end
PopSize = 50;
GenLimit = 200;
DisplayFlag = false;
MinCostBin40 = zeros(GenLimit+1, nMonte);
MinCostBin80 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostBin40(:,i) = ASCont(@Ackley, DisplayFlag, 40, PopSize, GenLimit); % # bins = 40, population size = 50
    MinCostBin80(:,i) = ASCont(@Ackley, DisplayFlag, 80, PopSize, GenLimit); % # bins = 80, population size = 50
end
MinCostBin40 = mean(MinCostBin40, 2);
MinCostBin80 = mean(MinCostBin80, 2);
SetPlotOptions
figure
plot(0:GenLimit, MinCostBin40, 'r--', 0:GenLimit, MinCostBin80, 'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('40 Intervals', '80 Intervals')