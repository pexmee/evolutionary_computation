function DEMonteDiff(nMonte)
% Compare DE performance using one or two difference vectors
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20; % number of Monte Carlo simulations
end
Problem = @Ackley;
Display = false;
GenLimit = 80;
LFlag = false; % Use the "/bin" option
Base = 1; % Use the best individual as the base vector
MinCostDiff1 = zeros(GenLimit+1, nMonte);
MinCostDiff2 = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostDiff1(:,i) = DE(Problem, Display, GenLimit, LFlag, Base, 1);
    MinCostDiff2(:,i) = DE(Problem, Display, GenLimit, LFlag, Base, 2);
end
MinCostDiff1 = mean(MinCostDiff1, 2);
MinCostDiff2 = mean(MinCostDiff2, 2);
SetPlotOptions
figure
plot(0:GenLimit,MinCostDiff2,'r--', 0:GenLimit,MinCostDiff1,'b-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('Two difference vectors', 'One difference vector')
