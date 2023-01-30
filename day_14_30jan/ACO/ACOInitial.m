function ACOInitial

NumAnts = 5000;
k = 20;
h = 2;
m1 = 0;
m2 = 0;
mArray = zeros(NumAnts+1, 2);
for i = 0 : NumAnts
    p1 = (m1 + k)^h / ((m1+k)^h + (m2+k)^h);
    if rand < p1
        m1 = m1 + 1;
    else
        m2 = m2 + 1;
    end
    mArray(i+1, :) = 100 * [m1 m2] / (i+1);
end
SetPlotOptions
figure
subplot(2, 1, 1), hold on
plot(0:100, mArray(1:101, 1), 'r-')
plot(0:100, mArray(1:101, 2), 'b--')
xlabel('Time (Number of Ants)')
ylabel('Percent of Ants')

subplot(2, 1, 2), hold on
plot(0:NumAnts, mArray(:, 1), 'r-')
plot(0:NumAnts, mArray(:, 2), 'b--')
legend('Path 1', 'Path 2')
xlabel('Time (Number of Ants)')
ylabel('Percent of Ants')

