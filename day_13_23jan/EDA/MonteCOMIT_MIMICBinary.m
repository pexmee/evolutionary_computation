function MonteCOMIT_MIMICBinary(ProblemFunction)
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @AckleyDisc;
end
nMonte = 20;
GenLimit = 20;
DisplayFlag = false;
MinCostMIMIC = zeros(GenLimit+1, nMonte);
MinCostCOMIT = zeros(GenLimit+1, nMonte);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    RandSeed = fix(sum(100*clock));
    MinCostMIMIC(:, i) = MIMICBinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, 40, false, false, true); 
    MinCostCOMIT(:, i) = MIMICBinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, 40, false, false, false); 
end
MinCostMIMIC = mean(MinCostMIMIC, 2);
MinCostCOMIT = mean(MinCostCOMIT, 2);
% Plot average of best cost at each generation
figure, hold on
plot(0:GenLimit, MinCostMIMIC, 'b-.')
plot(0:GenLimit, MinCostCOMIT, 'r-')
xlabel('Generation')
ylabel('Minimum Cost')
legend('MIMIC', 'COMIT')