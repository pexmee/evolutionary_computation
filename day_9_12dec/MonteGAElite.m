function MonteGAElite(nMonte)
% Explore the effect of elitism on a GA
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20;
end
GenLimit = 40;
Gray = false;
DisplayFlag = false;
PopSize = 20;
Pmutate = 0.02;
MinDomain = -5;
MaxDomain = +5;
Dimension = 2;
MinCostNoElite = zeros(GenLimit+1, nMonte);
MinCostYesElite = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    MinCostNoElite(:, i) = GA(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, Gray, 0, PopSize, ...
        Pmutate, MinDomain, MaxDomain, Dimension); 
    MinCostYesElite(:, i) = GA(@AckleyDisc, DisplayFlag, RandSeed, GenLimit, Gray, 2, PopSize, ...
        Pmutate, MinDomain, MaxDomain, Dimension); 
end
MinCostNoElite = mean(MinCostNoElite, 2);
MinCostYesElite = mean(MinCostYesElite, 2);
figure, hold on, box on
plot(0:GenLimit, MinCostNoElite, 'b-')
plot(0:GenLimit, MinCostYesElite, 'r--')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Binary Coding without Elitism', 'Binary Coding with Elitism')