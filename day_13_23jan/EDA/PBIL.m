function [MinCost] = PBIL(ProblemFunction, DisplayFlag, GenLimit, LearningRate, UpdateCount, SigmaFactor)

% Probability Based Incremental Learning (PBIL) for optimizing a continuous function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         LearningRate = ... well, learning rate - what did you expect?
%         UpdateCount = number of best and worst solutions used to update probability vector
%         SigmaFactor = initial and final proportion of std dev of solution generation to parameter range
%                       (two elements - one for initial, and one for final)
% OUTPUT: MinCost = array of best solution, one element for each generation

if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 100;
end
OPTIONS.Maxgen = GenLimit;
if ~exist('LearningRate', 'var') || isempty(LearningRate)
    LearningRate = 0.1;
end
if exist('UpdateCount', 'var') && ~isempty(UpdateCount)
    UpdateFromBest = UpdateCount; % # good individuals used to update the probability vector each generation
    UpdateFromWorst = UpdateCount; % # bad individuals used to update the probability vector each generation
else
    UpdateFromBest = 5; 
    UpdateFromWorst = 5;
end
if ~exist('SigmaFactor', 'var') || isempty(SigmaFactor)
    SigmaFactor = [1/10, 1/50];
end

OPTIONS.popsize = 50;
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
Range = OPTIONS.MaxDomain - OPTIONS.MinDomain; % parameter range

sigmaMax = SigmaFactor(1) * Range; % initial std dev of generated chromosomes
sigmaMin = SigmaFactor(2) * Range; % final std dev of generated chromosomes
shiftMutate = 0; % probability vector mutation shift magnitude
ProbVec = OPTIONS.MinDomain + 2 * Range .* rand(1, OPTIONS.numVar); % initial probability vector

% Begin the evolution loop
for GenIndex = 0 : OPTIONS.Maxgen
    % Determine std dev of generated chromosomes
    sigma = sigmaMax - GenIndex * (sigmaMax - sigmaMin) / OPTIONS.Maxgen;
    % Generate a population based on the probability vector.
    for popindex = OPTIONS.Keep+1 : OPTIONS.popsize
        chrom = ProbVec + sigma .* randn(1, OPTIONS.numVar);
        Population(popindex).chrom = max(OPTIONS.MinDomain, min(OPTIONS.MaxDomain, chrom));
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    % Sort from best to worst
    Population = PopSort(Population);
    % Probability vector update from best population members
    for k = 1 : UpdateFromBest
        ProbVec = ProbVec + LearningRate * (Population(k).chrom - ProbVec);
    end
    % Probability vector update from worst population members
    for k = OPTIONS.popsize-UpdateFromWorst+1 : OPTIONS.popsize
        ProbVec = ProbVec - LearningRate * (Population(k).chrom - ProbVec);
    end
    % Mutation of the probability vector
    for i = 1 : OPTIONS.numVar
        if rand < OPTIONS.pmutate
            ProbVec(i) = ProbVec(i) + shiftMutate * 2 * (rand - 0.5);
        end
    end
    ProbVec = max(OPTIONS.MinDomain, min(OPTIONS.MaxDomain, ProbVec));
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return