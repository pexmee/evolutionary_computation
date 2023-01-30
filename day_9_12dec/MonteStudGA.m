function [MeanMin, MeanMinNorm, BestMin, BestMinNorm] = MonteStudGA(nMonte)
% Monte Carlo execution of evolutionary algorithm software to explore the effect of the stud option
% OUTPUT MeanMin is the mean of the best solution found. It is an
%        nFunction x nBench array, where nFunction is the number of optimization
%        functions that are used, and nBench is the number of benchmarks that
%        are optimized.
% OUTPUT MeanMinNorm is MeanMin normalized to a minimum of 1 for each benchmark.
% OUTPUT BestMin is the best solution found by each optimization function
%        for each benchmark.
% OUTPUT BestMinNorm is BestMin normalized to a minimum of 1 for each benchmark.
if ~exist('nMonte', 'var') || isempty(nMonte)
    nMonte = 50; % number of Monte Carlo runs
end
% Benchmark functions
 Bench = [     %     multimodal? separable?  regular?
 'Ackley     '; %     y           n           y
 'Fletcher   '; %     y           n           n
 'Griewank   '; %     y           n           y
 'Penalty1   '; %     y           n           y
 'Penalty2   '; %     y           n           y
 'Quartic    '; %     n           y           y
 'Rastrigin  '; %     y           y           y
 'Rosenbrock '; %     n           n           y
 'Schwefel12 '; %     y           y           n
 'Schwefel221'; %     n           n           y
 'Schwefel222'; %     y           n           n
 'Schwefel226'; %     n           n           n
 'Sphere     '; %     n           y           y
 'Step       ']; %    n           y           n
nBench = size(Bench, 1);
MeanMin = zeros(2, nBench);
BestMin = inf(2, nBench);
for j = 1 : nBench
    disp(['Benchmark function ', num2str(j), '/', num2str(nBench)]);
    for k = 1 : nMonte
        % GA without the stud option
        [Cost] = eval(['GAStud(@', Bench(j,:), ', false, false);']);
        MeanMin(1,j) = ((k - 1) * MeanMin(1,j) + Cost(end)) / k;
        BestMin(1,j) = min(BestMin(1,j), Cost(end));
        % GA with the stud option
        [Cost] = eval(['GAStud(@', Bench(j,:), ', false, true);']);
        MeanMin(2,j) = ((k - 1) * MeanMin(2,j) + Cost(end)) / k;
        BestMin(2,j) = min(BestMin(2,j), Cost(end));
    end
end
% Normalize the results
if min(MeanMin) == 0
    MeanMinNorm = [];
else
    MeanMinNorm = MeanMin * diag(1./min(MeanMin));
end
if min(BestMin) == 0
    BestMinNorm = [];
else
    BestMinNorm = BestMin * diag(1./min(BestMin));
end