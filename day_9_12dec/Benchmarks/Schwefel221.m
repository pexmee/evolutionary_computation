function [InitFunction, CostFunction] = Schwefel221
InitFunction = @Init;
CostFunction = @Cost;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
OPTIONS.MinDomain = -100 * ones(1, OPTIONS.numVar);
OPTIONS.MaxDomain = +100 * ones(1, OPTIONS.numVar);
OPTIONS.OrderDependent = false; % the order of the independent variables does not matter
OPTIONS = ComputeRandomShift(OPTIONS);
% Initialize population
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    chrom = OPTIONS.MinDomain + (OPTIONS.MaxDomain - OPTIONS.MinDomain) .* rand(1,OPTIONS.numVar);
    Population(popindex).chrom = chrom;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, OPTIONS)
% Compute the cost of each member in Population
p = length(Population(1).chrom);
for popindex = 1 : length(Population)
    Population(popindex).cost = 0;
    for i = 1 : p
        x = Population(popindex).chrom(i) - OPTIONS.ShiftAmount(i);
        Population(popindex).cost = max(Population(popindex).cost, abs(x));
    end
end
return