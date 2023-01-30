function [Population] = ReplaceDupsTSP(Population, MutationMethod, DistanceArray)
% Replace duplicate individuals in the TSP population by mutating them.
% This routine does not make 100% certain that no duplicates exist in the population, but
% it mutates duplicate individuals, and then checks again and replaces duplicate individuals
% with random tours, so after this routine exits, there is only a small probability that 
% any duplicate individuals are in the population.
N = length(Population);
NumCities = length(Population(1).Tour) - 1;
% Mutate duplicate tours
for i = 1 : N
    Tour1 = Population(i).Tour;
    for j = i+1 : N
        Tour2 = Population(j).Tour;
        if isequal(Tour1, Tour2)
            Population(j).Tour = MutateTSP(Population(j).Tour, MutationMethod);
            Population(j).Distance = CalcDistance(Population(j).Tour, DistanceArray);
        end
    end
end
% Check again for duplicates, this time replacing duplicate tours with random tours
for i = 1 : N
    Tour1 = Population(i).Tour;
    for j = i+1 : N
        Tour2 = Population(j).Tour;
        if isequal(Tour1, Tour2)
            Population(j).Tour = randperm(NumCities);
            % Make the tour a round trip that begins and ends at city #1
            ndx = find(Population(j).Tour == 1, 1, 'first');
            Population(j).Tour = circshift(Population(j).Tour, [0, 1-ndx]);
            Population(j).Tour(end+1) = Population(j).Tour(1);
            Population(j).Distance = CalcDistance(Population(j).Tour, DistanceArray);
        end
    end
end
return