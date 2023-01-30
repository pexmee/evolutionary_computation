function [MinDistance] = AS(Filename, DisplayFlag)

% Ant colony optimization algorithm for solving a traveling salesman problem.

% INPUTS: Filename is the name of the file that has the TSP data (no extension; .TSP assumed).
%            Assumed format of file is from http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95.
%         DisplayFlag says whether or not to display information during iterations and plot results.
% OUTPUT: MinDistance = array of minimum distances obtained at each generation

if ~exist('Filename', 'var')
    Filename = 'ulysses16';
end
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
disp('Initializing ...');

% ACO parameter initialization
MaxGen = 20; % Generation limit
Keep = 2; % Elitism parameter: how many of the best individuals to keep from one generation to the next
tau0 = 1e-6; % initial pheromone value, between 0 and 0.5
Q = 20; % pheromonone update constant, between 0 and 100
rhog = 0.9; % global pheromone decay rate, between 0 and 1
alpha = 1; % pheromone sensitivity, between 1 and 5
beta = 5; % heuristic sensitivity, between 0 and 15

RandSeed = sum(100*clock);
if DisplayFlag
    disp(['Random seed = ', num2str(RandSeed)]);
end
rng(RandSeed);

% Read the input file and fill up the Coordinates array with the coordinates of each city
[Coord, EdgeWeightType] = GetCoordinates(Filename);
if isequal(EdgeWeightType, 'EUC_2D')
elseif isequal(EdgeWeightType, 'GEO')
    Coord = GetLongLat(Coord);
else
    disp('Edge Weight Type not recognized')
    return
end
DistanceArray = CreateDistanceArray(Coord, EdgeWeightType);

% Initialize the population and calculate the distance of each tour
NumCities = size(Coord, 1);
PopSize = NumCities + 1;
Population = struct('Tour', cell([1 PopSize]), 'Distance', cell([1 PopSize]));
for i = 1 : PopSize
    Population(i).Tour = randperm(NumCities);
    Population(i).Tour(end+1) = Population(i).Tour(1); % Make the tour a round trip
    Population(i).Distance = CalcDistance(Population(i).Tour, DistanceArray);
end
[Population] = PopSortTSP(Population);

MinDistance = zeros(1, MaxGen);
AvgDistance = zeros(1, MaxGen);
MinDistance(1) = Population(1).Distance;
AvgDistance(1) = mean([Population.Distance]);
if DisplayFlag
    disp(['The best and mean of Generation # 0 are ', ...
        num2str(MinDistance(1)), ' and ', num2str(AvgDistance(1))]);
    % Plot the initial best tour
    SetPlotOptions
    PlotTour(Coord, Population(1).Tour)
    title('Best Initial Tour')
end

tau = tau0 * ones(NumCities, NumCities); % initial pheromone values
for i = 1 : NumCities
    tau(i, i) = 0;
end
% Begin the optimization loop
for GenIndex = 1 : MaxGen
    % Use the probabilities to generate new solutions
    for k = Keep+1 : PopSize % for each non-elite ant
        for q = 1 : NumCities-1 % loop until all cities are visited
            City = Population(k).Tour(q); % Current city of the k-th ant
            % Calculate the denominator of the probability calculation
            Denom = 0;
            for  m = 1 : NumCities % for each potential next city
                if ~isempty(find(Population(k).Tour(1:q) == m, 1)), continue, end
                Denom = Denom + tau(City, m)^alpha / DistanceArray(City, m)^beta;
            end
            % Calculate the probability of travel to each potential city
            Prob = zeros(NumCities, 1); 
            for  j = 1 : NumCities % for each potential next city
                if ~isempty(find(Population(k).Tour(1:q) == j, 1)), continue, end
                Prob(j) = tau(City, j)^alpha / DistanceArray(City, j)^beta / Denom;
            end
            % Use the probabilities to decide which city to travel to next
            Select_index = 1;
            RandomNumber = rand;
            SelectProb = Prob(1);
            while SelectProb < RandomNumber
                Select_index = Select_index + 1;
                if Select_index >= NumCities, break, end
                SelectProb = SelectProb + Prob(Select_index);
            end
            Population(k).Tour(q+1) = Select_index;
        end
        Population(k).Tour(NumCities+1) = Population(k).Tour(1); % Make the tour a round trip
        Population(k).Distance = CalcDistance(Population(k).Tour, DistanceArray);
    end
    [Population] = PopSortTSP(Population);
    tau = (1 - rhog) * tau; % pheromone decay
    % Use each ant to update the pheromone between each pair of cities
    dtau = zeros(NumCities, NumCities); % number of ants that went from city i to city j
    for k = 1 : PopSize
        for i = 1 : NumCities
            FromCity = Population(k).Tour(i);
            ToCity = Population(k).Tour(i+1);
            dtau(FromCity, ToCity) = dtau(FromCity, ToCity) + Q / Population(k).Distance;
        end
    end
    tau = tau + dtau;
    MinDistance(GenIndex+1) = Population(1).Distance;
    AvgDistance(GenIndex+1) = mean([Population.Distance]);
    if DisplayFlag
        disp(['The best and mean of Generation # ', num2str(GenIndex), ' are ', ...
            num2str(MinDistance(GenIndex+1)), ' and ', num2str(AvgDistance(GenIndex+1))]);
    end
end
if DisplayFlag
    ConcludeTSP(Population, MinDistance, AvgDistance, Coord);
    PlotBestTour(Filename, Coord, DistanceArray);
end
return