function MonteESmulambdaAdaptAll
% Save Matlab figure files output from MonteESmulambdaAdapt.m for all benchmarks
MonteESmulambdaAdapt(@Ackley);
saveas(gcf, 'ESCommaAckley.fig');
MonteESmulambdaAdapt(@Fletcher);
saveas(gcf, 'ESCommaFletcher.fig');
MonteESmulambdaAdapt(@Griewank);
saveas(gcf, 'ESCommaGriewank.fig');
MonteESmulambdaAdapt(@Penalty1);
saveas(gcf, 'ESCommaPenalty1.fig');
MonteESmulambdaAdapt(@Penalty2);
saveas(gcf, 'ESCommaPenalty2.fig');
MonteESmulambdaAdapt(@Quartic);
saveas(gcf, 'ESCommaQuartic.fig');
MonteESmulambdaAdapt(@Rastrigin);
saveas(gcf, 'ESCommaRastrigin.fig');
MonteESmulambdaAdapt(@Rosenbrock);
saveas(gcf, 'ESCommaRosenbrock.fig');
MonteESmulambdaAdapt(@Schwefel12);
saveas(gcf, 'ESCommaSchwefel12.fig');
MonteESmulambdaAdapt(@Schwefel221);
saveas(gcf, 'ESCommaSchwefel221.fig');
MonteESmulambdaAdapt(@Schwefel222);
saveas(gcf, 'ESCommaSchwefel222.fig');
MonteESmulambdaAdapt(@Schwefel226);
saveas(gcf, 'ESCommaSchwefel226.fig');
MonteESmulambdaAdapt(@Sphere);
saveas(gcf, 'ESCommaSphere.fig');
MonteESmulambdaAdapt(@Step);
saveas(gcf, 'ESCommaStep.fig');