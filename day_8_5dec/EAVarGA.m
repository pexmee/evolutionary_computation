function [MinCost, AvgCost] = EAVarGA(ProblemFunction, DisplayFlag, RandSeed, GenLimit, GrayFlag, NumElites)

% Genetic algorithm for optimizing a general function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%            the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = whether or not to display information during iterations and plot results.
%         RandSeed = random number seed
%         GenLimit = number of generations
%         GrayFlag = whether to use gray coding or binary coding
%         NumElites = number of elite individuals each generation

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @WorstCaseProblem;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('RandSeed', 'var') || isempty(RandSeed)
    RandSeed = round(sum(100*clock));
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 50;
end
OPTIONS.Maxgen = GenLimit;
if ~exist('GrayFlag', 'var') || isempty(GrayFlag)
    GrayFlag = true;
end
OPTIONS.Gray = GrayFlag;
if ~exist('NumElites', 'var') || isempty(NumElites)
    NumElites = 2;
end
OPTIONS.Keep = NumElites;
OPTIONS.clearDups = false;

[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed);

Xover_Type = 1; % crossover type: 1 = single point, 2 = two point, 3 = uniform
pcross = 1; % crossover probability
child = zeros(OPTIONS.popsize-OPTIONS.Keep, OPTIONS.numVar);

% Begin the evolution loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Compute fitness, which increases as cost decreases.
    Fitness = max([Population.cost]) + 1 - [Population.cost];
    FitnessSum = sum(Fitness);
    for k = 1 : 2 : OPTIONS.popsize-OPTIONS.Keep % begin selection/crossover loop
        % Select two parents to mate and create two children - roulette wheel selection
        mate = [0 0];
        for selParents = 1 : 2
            Random_Cost = rand * FitnessSum;
            Select_Cost = Fitness(1);
            Select_index = 1;
            while Select_Cost < Random_Cost 
                Select_index = Select_index + 1;
                if Select_index >= OPTIONS.popsize, break, end
                Select_Cost = Select_Cost + Fitness(Select_index);
            end
            mate(selParents) = Select_index;
        end
        Parent(1, :) = Population(mate(1)).chrom;
        Parent(2, :) = Population(mate(2)).chrom;
        % Crossover
        switch Xover_Type
            case 1
                % single point crossover
                if rand < pcross
                    % crossover the parents
                    Xover_Pt = ceil(rand * OPTIONS.numVar);
                    child(k, :) = [Parent(1, 1:Xover_Pt), Parent(2, Xover_Pt+1:OPTIONS.numVar)];
                    child(k+1, :) = [Parent(2, 1:Xover_Pt), Parent(1, Xover_Pt+1:OPTIONS.numVar)];
                else
                    child(k, :) = Parent(1, :);
                    child(k+1, :) = Parent(2, :);
                end
            case 2
                % two point crossover
                if rand < pcross
                    Xover_Pt1 = ceil(rand * OPTIONS.numVar);
                    Xover_Pt2 = ceil(rand * OPTIONS.numVar);
                    if Xover_Pt1 > Xover_Pt2
                        temp = Xover_Pt2;
                        Xover_Pt2 = Xover_Pt1;
                        Xover_Pt1 = temp;
                    end
                    child(k, :) = [Parent(1, 1:Xover_Pt1) Parent(2, Xover_Pt1+1:Xover_Pt2) Parent(1, Xover_Pt2+1:OPTIONS.numVar)];
                    child(k+1, :) = [Parent(2, 1:Xover_Pt1) Parent(1, Xover_Pt1+1:Xover_Pt2) Parent(2, Xover_Pt2+1:OPTIONS.numVar)];
                else
                    child(k, :) = Parent(1, :);
                    child(k+1, :) = Parent(2, :);
                end
            case 3 
                % uniform crossover
                if rand < pcross
                    for i = 1 : OPTIONS.numVar
                        if rand < 0.5
                            child(k, i) = Parent(1, i);
                            child(k+1, i) = Parent(2, i);
                        else
                            child(k, i) = Parent(2, i);
                            child(k+1, i) = Parent(1, i);
                        end
                    end
                else
                    child(k, :) = Parent(1, :);
                    child(k+1, :) = Parent(2, :);
                end
        end
    end % end selection/crossover loop
    % Replace the non-elite population members with the new children
    for k = OPTIONS.Keep+1 : 2 : OPTIONS.popsize
        Population(k).chrom = child(k-OPTIONS.Keep, :);
        Population(k+1).chrom = child(k-OPTIONS.Keep+1, :);
    end
    % Mutation - but don't allow the elites to be mutated
    for individual = OPTIONS.Keep+1 : OPTIONS.popsize 
        for parnum = 1 : OPTIONS.numVar
            if OPTIONS.pmutate > rand
                Population(individual).chrom(parnum) = ~Population(individual).chrom(parnum);
            end
        end
    end
    % Make sure the population does not have duplicates.
    if OPTIONS.clearDups
        Population = ClearDups(Population);
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