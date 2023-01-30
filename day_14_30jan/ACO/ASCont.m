function [MinCost] = ASCont(ProblemFunction, DisplayFlag, NumBins, PopSize, GenLimit, ...
    tauMin, tauMax, NumBest, phi, q0)

% Ant system algorithm for optimizing a general continuous function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%            the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = whether or not to display information to the screen and plot results.
%         NumBins = number of bins per domain
%         PopSize = population size
%         GenLimit = generation limit
%         tauMin, tauMax = minimum and maximum pheromone for MMAS
%         NumBest = how many ants are allowed to deposit pheromone each generation
%         phi = local pheromone decay rate, between 0 (standard AS) and 1
%         q0 = ACS exploration constant, between 0 (standard AS) and 1
% OUTPUT: MinCost = array of best cost values at each generation

OPTIONS = [];
if ~exist('ProblemFunction', 'var')
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var')
    DisplayFlag = true;
end
if ~exist('NumBins', 'var')
    NumBins = 100; % discretization levels of search domain for pheromone deposits
end
if ~exist('PopSize', 'var')
    PopSize = 40;
end
if ~exist('GenLimit', 'var')
    GenLimit = 100;
end
if ~exist('tauMin', 'var')
    tauMin = 0;
end
if ~exist('tauMax', 'var')
    tauMax = inf;
end
if ~exist('NumBest', 'var')
    NumBest = inf;
end
if ~exist('phi', 'var')
    phi = 0;
end
if ~exist('q0', 'var')
    q0 = 0;
end
OPTIONS.popsize = PopSize;
OPTIONS.Maxgen = GenLimit;
OPTIONS.numVar = 20; % problem dimension

[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);

% ACO parameter initialization
tau0 = 1e-6; % initial pheromone value, between 0 and 0.5
Q = 20; % pheromonone update constant, between 0 and 100
rhog = 0.9; % global pheromone decay rate, between 0 and 1
alpha = 1; % pheromone sensitivity, between 1 and 5
tau = cell(OPTIONS.numVar, 1); % cell array of pheromone values, one cell array element per problem dimension
for i = 1 : OPTIONS.numVar
    tau{i}(1, :) = linspace(OPTIONS.MinDomain(i), OPTIONS.MaxDomain(i), NumBins+1);
    tau{i}(2, :) = tau0 * ones(1, NumBins+1);
    tau{i}(2, end) = 0;
end
DeltaTau = tau{1}(1,2) - tau{1}(1,1);

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    for i = 1 : OPTIONS.numVar
        tau{i}(2,:) = (1 - rhog) * tau{i}(2,:); % Global pheromone decay
    end
    % Use each solution to update the pheromone for each parameter value
    for k = 1 : min(NumBest, OPTIONS.popsize)
        Cost = Population(k).cost;
        Chrom = Population(k).chrom;
        for i = 1 : length(Chrom)
            Bin = find(tau{i}(1,:) > Chrom(i), 1) - 1;
            Bin = min(max(Bin, 1), NumBins);
            if Cost <= 0
                tau{i}(2,Bin) = max(tau{i}(2,:));
            else
                tau{i}(2,Bin) = tau{i}(2,Bin) + Q / Cost;
            end
        end    
    end
    for i = 1 : OPTIONS.numVar
        tau{i}(2,:) = min( max( tau{i}(2,:), tauMin), tauMax);
    end
    % Use pheromone-based probabilities to generate new solutions
    for k = OPTIONS.Keep+1 : OPTIONS.popsize
        for j = 1 : OPTIONS.numVar
            if OPTIONS.pmutate > rand
                % Mutation
                Population(k).chrom(j) = OPTIONS.MinDomain(j) + ...
                    (OPTIONS.MaxDomain(j) - OPTIONS.MinDomain(j)) * rand;
            else
                % Generate probabilities based on pheromone amounts
                Prob = tau{j}(2,:) .^ alpha;
                Prob = Prob / sum(Prob);
                [~, Select_index] = max(Prob);
                if rand > q0
                    SelectProb = Prob(1);
                    Select_index = 1;
                    RandomNumber = rand;
                    while SelectProb < RandomNumber
                        Select_index = Select_index + 1;
                        if Select_index >= NumBins, break, end
                        SelectProb = SelectProb + Prob(Select_index);
                    end
                end
                Population(k).chrom(j) = tau{j}(1,Select_index) + DeltaTau * rand;
                % local pheromone update
                tau{j}(2,Select_index) = (1 - phi) * tau{j}(2,Select_index) + phi * tau0;
            end
       end
    end
    if OPTIONS.clearDups
        % Make sure the population does not have duplicates. 
        Population = ClearDups(Population, OPTIONS);
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Sort from best to worst
    Population = PopSort(Population);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return
