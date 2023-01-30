function ACSContMonte1(nMonte)
% Monte Carlo ant system simulation to explore the effect of the local pheromone decay constant
% INPUT: nMonte = number of Monte Carlo simulations
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20;
end
NumIntervals = 20;
PopSize = 40;
GenLimit = 100;
tauMin = 0;
tauMax = inf;
NumBest = 4;
q0 = 0;
DisplayFlag = false;
MinCostphi0pt000 = zeros(GenLimit+1, nMonte);
MinCostphi0pt001 = zeros(GenLimit+1, nMonte);
MinCostphi0pt010 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostphi0pt000(:,i) = ASCont(@Ackley, DisplayFlag, NumIntervals, PopSize, GenLimit, ...
        tauMin, tauMax, NumBest, 0.000, q0);
    MinCostphi0pt001(:,i) = ASCont(@Ackley, DisplayFlag, NumIntervals, PopSize, GenLimit, ...
        tauMin, tauMax, NumBest, 0.001, q0);
    MinCostphi0pt010(:,i) = ASCont(@Ackley, DisplayFlag, NumIntervals, PopSize, GenLimit, ...
        tauMin, tauMax, NumBest, 0.010, q0);
end
MinCostphi0pt000 = mean(MinCostphi0pt000, 2);
MinCostphi0pt001 = mean(MinCostphi0pt001, 2);
MinCostphi0pt010 = mean(MinCostphi0pt010, 2);
SetPlotOptions
figure
plot(0:GenLimit, MinCostphi0pt000, 'r--', 0:GenLimit, MinCostphi0pt001, 'b-', ...
    0:GenLimit, MinCostphi0pt010, 'k:')
xlabel('Generation')
ylabel('Minimum Cost')
legend('\phi = 0', '\phi = 0.001', '\phi = 0.01')