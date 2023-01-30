function DEMonteF(nMonte, Problem)
% Compare DE performance using dithered, jittered, or constant F
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 100; % number of Monte Carlo simulations
end
if ~exist('Problem', 'var') || isempty(Problem)
    Problem = @Fletcher;
end
Display = false;
GenLimit = 30;
LFlag = false; % Use the "/bin" option
Base = 1; % Use the best individual as the base vector
NumDiff = 1; % use 1 difference vector
MinCostF0 = zeros(GenLimit+1, nMonte);
MinCostF1 = zeros(GenLimit+1, nMonte);
MinCostF2 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostF0(:,i) = DE(Problem, Display, GenLimit, LFlag, Base, NumDiff, 0);
    MinCostF1(:,i) = DE(Problem, Display, GenLimit, LFlag, Base, NumDiff, 1);
    MinCostF2(:,i) = DE(Problem, Display, GenLimit, LFlag, Base, NumDiff, 2);
end
MinCostF0 = mean(MinCostF0, 2);
MinCostF1 = mean(MinCostF1, 2);
MinCostF2 = mean(MinCostF2, 2);
SetPlotOptions
figure
plot(0:GenLimit,MinCostF1,'b-', 0:GenLimit,MinCostF0,'r--', 0:GenLimit,MinCostF2,'k:')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Dithered F', 'Constant F', 'Jittered F')
