function DEMonteLbin(nMonte)
% Compare the "/L" and the "/bin" versions of DE
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20; % number of Monte Carlo simulations
end
Problem = @Ackley;
Display = false;
GenLimit = 80;
Basis = 1; % use the best individual as the basis vector
NumDiff = 1; % use 1 difference vector
MinCostLfalse = zeros(GenLimit+1, nMonte);
MinCostLtrue = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    MinCostLtrue(:,i) = DE(Problem, Display, GenLimit, true, Basis, NumDiff);
    MinCostLfalse(:,i) = DE(Problem, Display, GenLimit, false, Basis, NumDiff);
end
MinCostLfalse = mean(MinCostLfalse, 2);
MinCostLtrue = mean(MinCostLtrue, 2);
SetPlotOptions
figure
plot(0:GenLimit, MinCostLtrue, 'b-', 0:GenLimit, MinCostLfalse, 'r--')
xlabel('Generation')
ylabel('Minimum Cost')
legend('DE/best/1/L', 'DE/best/1/bin')
