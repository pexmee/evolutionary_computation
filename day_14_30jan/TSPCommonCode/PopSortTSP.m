function [Population, indices] = PopSortTSP(Population)
% Sort the TSP population members from best to worst
popsize = length(Population);
[Distances, indices] = sort([Population.Distance], 'ascend');
Tours = zeros(popsize, length(Population(1).Tour));
for i = 1 : popsize
    Tours(i, :) = Population(indices(i)).Tour;
end
for i = 1 : popsize
    Population(i).Tour = Tours(i, :);
    Population(i).Distance = Distances(i);
end
return