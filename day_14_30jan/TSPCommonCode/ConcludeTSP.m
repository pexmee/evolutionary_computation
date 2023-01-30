function ConcludeTSP(Population, MinDistance, AvgDistance, Coord)
% Finish up the TSP run by displaying and plotting some data.
% INPUTS: Population = final population of TSP candidate solutions
%         MinDistance = array of best solutions, one element per generation
%         AvgDistance = array of average solutions, one element per generation
%         Coord = n x 2 array of lat/long coordinates
Tour = Population(1).Tour; % best individual
disp(['Best tour found: ', num2str(Tour)]); 
figure
plot((0:length(AvgDistance)-1), AvgDistance, 'r--', (0:length(AvgDistance)-1), MinDistance, 'b-')
legend('Average Distance', 'Shortest Distance')
xlabel('Generation')
ylabel('Distance')
PlotTour(Coord, Tour), title('Best Tour Found')
% Check how many duplicate individuals there are in the population
N = length(Population);
DupFlag = zeros(1, N);
for i = 1 : N
    if DupFlag(i), continue, end
    Tour1 = Population(i).Tour;
    for j = i+1 : N
        Tour2 = Population(j).Tour;
        if isequal(Tour1, Tour2)
            DupFlag(i) = 1;
            DupFlag(j) = 1;
        end
    end
end
NumUnique = length(find(DupFlag == 0));
disp([num2str(NumUnique), ' unique individuals in final population.']);
return