function ACSContMonte2(nMonte)
% Monte Carlo ant system simulation to explore the effect of the exploration constant
% INPUT: nMonte = number of Monte Carlo simulations
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 100;
end
NumIntervals = 20;
PopSize = 40;
GenLimit = 50;
tauMin = 0;
tauMax = inf;
NumBest = 4;
phi = 0;
DisplayFlag = false;
q0Arr = [0, 0.001, 0.01, 0.1, 0.2, 0.5, 0.9];
numq0 = length(q0Arr);
MinCost = zeros(GenLimit+1, nMonte, numq0);
for i = 1 : nMonte
    disp(['Run # ', num2str(i), ' of ', num2str(nMonte)]);
    for q = 1 : numq0
        MinCost(:,i,q) = ASCont(@Ackley, DisplayFlag, NumIntervals, PopSize, GenLimit, ...
        tauMin, tauMax, NumBest, phi, q0Arr(q));
    end
end
MinCost = mean(MinCost, 2);
colors = ['k: '; 'r- '; 'b--'; 'm-.'; 'ko '; 'ro '; 'bo '];
ncolors = size(colors, 1);
SetPlotOptions
figure, hold on
Lgnd = cell(numq0, 1);
for q = 1 : numq0
    colorndx = mod(q, ncolors);
    if colorndx == 0
        colorndx = ncolors;
    end
    plot(0:GenLimit, MinCost(:,q), colors(colorndx,:))
    Lgnd{q} = ['q0 = ', num2str(q0Arr(q))];
end
legend(Lgnd);
xlabel('Generation')
ylabel('Minimum Cost')