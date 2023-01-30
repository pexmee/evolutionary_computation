function [MinCost] = UMDABinary(ProblemFunction, DisplayFlag, RandSeed, GenLimit, UpdatePct)

% Univariate marginal distribution algorithm (UMDA) for minimizing a function on a binary domain.
% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         RandSeed = random number seed for reproducibility of results
%         GenLimit = generation count
%         UpdatePct = percent of best solutions used to update probability vector
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
if ~exist('UpdatePct', 'var') || isempty(UpdatePct)
    UpdatePct = 10;
end
UpdatePct = UpdatePct / 100;

OPTIONS.Dim = 20;

OPTIONS.Dim = 5;
OPTIONS.Maxgen = 20;

OPTIONS.popsize = 100;
OPTIONS.pmutate = 0;
OPTIONS.Gray = false; 
OPTIONS.clearDups = false;

[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS, RandSeed);

% Begin the evolution loop
for GenIndex = 1 : OPTIONS.Maxgen
    ElitePop = Population(1 : OPTIONS.Keep);
    % Collect statistics
    BitProb = zeros(1, OPTIONS.numVar);
    Mcount = round(OPTIONS.popsize*UpdatePct);
    for i = 1 : Mcount
        BitProb = BitProb + Population(i).chrom;
    end
    BitProb = BitProb / Mcount;            
    % Generate a population based on the probability vector.
    for popindex = OPTIONS.Keep+1 : OPTIONS.popsize
        Population(popindex).chrom = (rand(1, OPTIONS.numVar) < BitProb);
    end
    % Mutation (don't mutate the elites)
    for i = OPTIONS.Keep+1 : OPTIONS.popsize
        for k = 1 : OPTIONS.numVar
            if rand < OPTIONS.pmutate
                Population(i).chrom(k) = ~Population(i).chrom(k);
            end
        end
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Sort from best to worst
    Population = PopSort(Population);
    % Replace the worst individuals with the best from the previous generation
    for i = 1 : OPTIONS.Keep
        Population(OPTIONS.popsize-i+1) = ElitePop(i);
    end
    Population = PopSort(Population);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return