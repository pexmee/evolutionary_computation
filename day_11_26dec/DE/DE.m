function [MinCost] = DE(ProblemFunction, DisplayFlag, GenLimit, LFlag, Base, NumDiff, FChange)

% Differential evolution algorithm for minimizing a continuous function.

% INPUTS: ProblemFunction = the handle of the function that returns 
%                           the handles of the initialization and cost functions.
%         DisplayFlag = true or false, whether or not to display and plot results.
%         GenLimit = generation limit
%         LFlag = whether to use the "/L" variation instead of the "/bin" variation
%         Base = whether to use best (1), random (2), or current (3) individual as the base vector
%         NumDiff = number of difference vectors (1 or 2)
%         FChange = whether to use constant F (0), dithered F (1), or jittered F (2)
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
if ~exist('LFlag', 'var') || isempty(LFlag)
    LFlag = false;
end
if ~exist('Base', 'var') || isempty(Base)
    Base = 1;
end
if ~exist('NumDiff', 'var') || isempty(NumDiff)
    NumDiff = 1;
end
if ~exist('FChange', 'var') || isempty(FChange)
    FChange = 1;
end

% Initialization
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
TrialVec = Population;

% DE parameter initialization
F0 = 0.4; % weighting factor
if NumDiff > 1
    F0 = F0 / 2;
end
F = F0;
CR = 0.9; % crossover constant

% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    for k = 1 : OPTIONS.popsize
        % Generate the mutant vector v
        if Base == 1
            r(1) = 1; % use the best individual as the base vector
        elseif Base == 2
            r(1) = randi(OPTIONS.popsize); % use a random individual as the base vector
        else
            r(1) = k;
        end
        while true
            r(2) = randi(OPTIONS.popsize);
            if (r(2) ~= r(1)), break, end
        end
        while true
            r(3) = randi(OPTIONS.popsize);
            if (r(3) ~= r(1)) && (r(3) ~= r(2)), break, end
        end
        if FChange == 1
            F = F0 + 0.4 * (rand - 0.5); % Dither F
        elseif FChange == 2
            F = F0 + 0.4 * (rand(1, OPTIONS.numVar) - 0.5); % Jitter F
        end
        vchrom = Population(r(1)).chrom + F .* (Population(r(2)).chrom - Population(r(3)).chrom);
        if NumDiff > 1
            while true
                r(4) = randi(OPTIONS.popsize);
                if (r(4) ~= r(1)) && (r(4) ~= r(2)) && (r(4) ~= r(3)), break, end
            end
            while true
                r(5) = randi(OPTIONS.popsize);
                if (r(5) ~= r(1)) && (r(5) ~= r(2)) && (r(5) ~= r(3)) && (r(5) ~= r(4)), break, end
            end
            vchrom = vchrom + F * (Population(r(4)).chrom - Population(r(5)).chrom);
        end
        % Generate the trial vector
        if LFlag % "/L" option
            L = randi(OPTIONS.numVar);
            s = randi(OPTIONS.numVar);
            TrialVec(k).chrom = Population(k).chrom;
            for j = s : min(s+L-1, OPTIONS.numVar)
                TrialVec(k).chrom(j) = vchrom(j);
            end
            if (s + L - 1) > OPTIONS.numVar
                for j = 1 : mod(s+L-1, OPTIONS.numVar)
                    TrialVec(k).chrom(j) = vchrom(j);
                end
            end 
        else % "/bin" option
            jrand = randi(OPTIONS.numVar);
            for j = 1 : OPTIONS.numVar
                if (rand < CR) || (j == jrand)
                    TrialVec(k).chrom(j) = vchrom(j);
                else
                    TrialVec(k).chrom(j) = Population(k).chrom(j);
                end
            end
        end
    end
    % Decide whether to keep the trial vector or original vector
    TrialVec = OPTIONS.CostFunction(TrialVec, OPTIONS);
    for k = 1 : OPTIONS.popsize
        if TrialVec(k).cost < Population(k).cost
            Population(k) = TrialVec(k);
        end
    end
    if OPTIONS.clearDups
        % Make sure the population does not have duplicates
        Population = ClearDups(Population, OPTIONS);
    end
    % Calculate cost
    Population = OPTIONS.CostFunction(Population, OPTIONS);
    Population = PopSort(Population);
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
return
