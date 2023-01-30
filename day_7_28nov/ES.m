function [MinCost] = ES(ProblemFunction, DisplayFlag, PlusFlag, AdaptFlag, Mu, Lambda, GenLimit, cAdapt)
% Evolutionary Strategy for optimizing a general continuous function.
%         Use AdaptFlag = 1 and GenLimit = 500 to reproduce Figure 6.4.
%         Use AdaptFlag = 2 to reproduce Figure 6.15.
% INPUTS: ProblemFunction = handle of the function that returns 
%         the handles of the initialization, cost, and feasibility functions.
%         DisplayFlag = whether or not to display information during iterations and plot results.
%         PlusFlag = whether to use the (mu+lambda) or (mu, lambda) ES
%         AdaptFlag = 0 for no adaptation of mutuation,
%                     1 to adapt the mutation rate using the 1/5 rule,
%                     2 to adapt the mutation rate randomly 
%         Mu = population size
%         Lambda = number of offspring
%         GenLimit = number of generations
%         cAdapt = adaptation constant c
% OUTPUT: MinCost = array of best cost, one element per generation
if ~exist('ProblemFunction', 'var') || isempty(ProblemFunction)
    ProblemFunction = @Ackley;
end
if ~exist('DisplayFlag', 'var') || isempty(DisplayFlag)
    DisplayFlag = true;
end
if ~exist('PlusFlag', 'var') || isempty(PlusFlag)
    PlusFlag = true;
end
if ~exist('AdaptFlag', 'var') || isempty(AdaptFlag)
    AdaptFlag = 0;
end
if ~exist('Mu', 'var') || isempty(Mu)
    Mu = 1;
end
if ~exist('Lambda', 'var') || isempty(Lambda)
    Lambda = 1;
end
if ~exist('GenLimit', 'var') || isempty(GenLimit)
    GenLimit = 100;
end
if ~exist('cAdapt', 'var') || isempty(cAdapt)
    cAdapt = 0.817;
end
Mutation.AdaptConst = cAdapt; % adaptation constant for mutation standard deviation
% ES parameters
OPTIONS.popsize = Mu;
OPTIONS.Maxgen = GenLimit; % generation count limit
OPTIONS.numVar = 20; % problem dimension
% Initialization
[OPTIONS, MinCost, AvgCost, Population, MinConstrViol, AvgConstrViol] = ...
    Init(DisplayFlag, ProblemFunction, OPTIONS);
if ~PlusFlag && (Lambda < OPTIONS.popsize)
    disp('ERROR - If using (mu, lambda)-ES, then lambda must be >= mu');
    return
end
if (AdaptFlag == 1) && (OPTIONS.popsize > 1)
    disp('ERROR - It does not make sense for AdaptFlag=1 if mu>1');
    return
end
% Initialize the standard deviations proportional to the standard deviation of the initial population,
% assuming that the initial population is uniformly distributed in its domain
StdInit = (OPTIONS.MaxDomain - OPTIONS.MinDomain) / 2 / sqrt(3);
for k = 1 : OPTIONS.popsize
    Population(k).std = StdInit / 100;
end
Children = Population(1 : min(Lambda, OPTIONS.popsize)); % initialize the Children structure
SigmaArray = zeros(OPTIONS.Maxgen+1, OPTIONS.numVar);
SigmaArray(1, :) = Population(1).std;
SuccessArray = zeros(1, OPTIONS.Maxgen+1);
SuccessRate = 0;
tau = 1 / sqrt(2 * sqrt(OPTIONS.numVar)); % factor for adaptive mutation
tauprime = 1 / sqrt(2 * OPTIONS.numVar); % factor for adaptive mutation
% ES initialization
Mutations.Total = 0; % total number of mutations
Mutations.Successful = 0; % number of successful mutations
Mutations.Limit = 20; % number of mutations after which to check successful mutation rate
% Begin the optimization loop
for GenIndex = 1 : OPTIONS.Maxgen
    % Produce children via uniform crossover
    for k = 1 : Lambda
        if OPTIONS.popsize > 1
            % Randomly select two parents
            p(1) = randi([1 OPTIONS.popsize]);
            p(2) = randi([1 OPTIONS.popsize]);
            for i = 1 : OPTIONS.numVar
                Children(k).chrom(i) = Population(p(randi([1 2]))).chrom(i);
                Children(k).std(i) = Population(p(randi([1 2]))).std(i);
            end
        else
            % Population size = 1 - copy the single parent to the child
            Children(k) = Population(1);
        end
        if AdaptFlag == 2
            Children(k).std = Children(k).std .* exp(tauprime * randn + tau * randn(1, OPTIONS.numVar));
        end
        Children(k).chrom = Children(k).chrom + Children(k).std .* randn(1, OPTIONS.numVar); % Mutation
    end
    % Make sure the child population does not contain duplicates. 
    if OPTIONS.clearDups
        Children = ClearDups(Children, OPTIONS);
    end
    Children = OPTIONS.CostFunction(Children, OPTIONS); % Calculate the cost of the Children
    if PlusFlag
        TotalPopulation = [Population Children]; % (mu+lambda)-ES
        if AdaptFlag == 1 % 1/5 rule - keep track of the # of mutations, and the # of successful mutations
            Mutations.Total = Mutations.Total + Lambda;
            Mutations.Successful = Mutations.Successful + sum([Children.cost] < Population(1).cost);
        end
    else
        TotalPopulation = Children; % (mu,lambda)-ES
    end
    TotalPopulation = PopSort(TotalPopulation);
    Population = TotalPopulation(1 : OPTIONS.popsize);
    if AdaptFlag == 1 % 1/5 rule
        if Mutations.Total >= Mutations.Limit
            SuccessRate = Mutations.Successful / Mutations.Total;
            if SuccessRate > 1/5
                Population(1).std = Population(1).std / Mutation.AdaptConst;
            elseif SuccessRate < 1/5
                Population(1).std = Population(1).std * Mutation.AdaptConst;
            end
            Mutations.Total = 0;
            Mutations.Successful = 0;
        end
        SuccessArray(GenIndex+1) = SuccessRate;
    end     
    SigmaArray(GenIndex+1, :) = Population(1).std;
    [MinCost, AvgCost, MinConstrViol, AvgConstrViol] = ComputeCostAndConstrViol(Population, ...
        MinCost, AvgCost, MinConstrViol, AvgConstrViol, GenIndex, DisplayFlag);
end 
Conclude(DisplayFlag, OPTIONS, Population, MinCost, AvgCost, MinConstrViol, AvgConstrViol);
if DisplayFlag
    if AdaptFlag == 1
        figure, hold on
        plot(0:OPTIONS.Maxgen, SigmaArray(:,1), 'b-', 0:OPTIONS.Maxgen, SuccessArray, 'r--')
        plot([0 OPTIONS.Maxgen], [0.2 0.2], 'g:')
        xlabel('Generation')
        legend('Standard Deviation', 'Success Rate')
    elseif AdaptFlag == 2
        figure
        bar(1:OPTIONS.numVar, SigmaArray(end,:), 0.2, 'c')
        xlabel('Feature Index')
        ylabel('Std Dev of Best Individual')
        axisTemp = axis;
        axis([0 OPTIONS.numVar+1 axisTemp(3) axisTemp(4)])
    end
end
return