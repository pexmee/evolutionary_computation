function [MinCost] = cGABinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, alpha, PopSize)

% Compact genetic algorithm (cGA) for minimizing a function on a binary domain.
% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         RandSeed = random number seed for reproducibility of results
%         GenLimit = generation count
%         alpha = probability vector update size
% OUTPUT: MinCost = array of best solutions, one element for each generation

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @AckleyDisc;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('RandSeed', 'var') || isempty(RandSeed)
    RandSeed = [];
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 100;
end
OPTIONS.Maxgen = GenLimit;
if ~exist('alpha', 'var') || isempty(alpha)
    alpha = 0.01;
end
if ~exist('PopSize', 'var') || isempty(PopSize)
    PopSize = 2;
end
OPTIONS.popsize = PopSize;

OPTIONS.Dim = 20;
OPTIONS.Keep = 1; % 0 or 1 for this algorithm
OPTIONS.pmutate = 0;
OPTIONS.Gray = false; 
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed);
ProbVector = 0.5 * ones(1, OPTIONS.numVar);
ProbMin = 0.05;
ProbMax = 1 - ProbMin;

% Begin the evolution loop
for GenIndex = 1 : OPTIONS.Maxgen
    ElitePop = Population(1);
    % Generate a population based on the probability vector.
    for popindex = 1 : OPTIONS.popsize
        Population(popindex).chrom = (rand(1, OPTIONS.numVar) < ProbVector);
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Sort from best to worst
    if OPTIONS.Keep > 0
        Population = PopSort([Population, ElitePop]);
        Population = Population(1 : OPTIONS.popsize);
    else
        Population = PopSort(Population);
    end
    % Update the probability vector
    for i = 1 : OPTIONS.numVar
        if Population(1).chrom(i) ~= Population(end).chrom(i)
            if Population(1).chrom(i) > 0
                ProbVector(i) = ProbVector(i) + alpha;
                ProbVector(i) = max(min(ProbVector(i), ProbMax), ProbMin);
            end
        end
    end
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
if ElitePop.cost < Population(1).cost
    Population(1) = ElitePop;
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
if DisplayFlag
    disp(['Probability vector = ', num2str(ProbVector)]); 
end
return