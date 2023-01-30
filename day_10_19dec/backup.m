function PSONeighborsMonte(nMonte, Genlimit)
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 20;
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 40;
end
c1 = 2;
c2 = 2;
c3 = 2;
K = 0.9;
NeighborsArr = [0, 5, 10];
numOptions = length(NeighborsArr);
MinCost = zeros(numOptions, GenLimit+1);
ColorLine = ['r--';'k: ';'b- ';'g-.';'m--'];
for i = 1: nMonte
    disp(['Run #', num2str(i), ' of ', num2str(nMonte)]);
    for r = 1 : numOptions
        MinCost(r, :, i) = PSO(@Sphere, false, c1,c2,c3,K,GenLimit,NeighborsArr(r));
    end
end
MinCost = mean(MinCost,3);
SetPlotOptions
figure, hold on
LegendArray = cell(numOptions,1);
for r = 1: numOptions
    plot(0:GenLimit, MinCost(r,:), ColorLine(r,:))
    LegendArray{r} = ['Neighborhood size = ', num2str(NeighborsArr(r))];
end
xlabel('Generation')
ylabel('Minimum Cost')
legend(LegendArray)
box on
