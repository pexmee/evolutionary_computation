function Kullback

% Reproduce the Kullback-Leibler example in the Chow 1968 paper.
% Note that entropy is calculated using the natural log.

% x corresponds to a matrix of candidate solutions to an optimization problem.
% One row for each individual, one column for each bit in the search domain.
% Individuals are distributed according to the first two columns of Table I.
x =[0 0 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 1;
    0 0 1 0;
    0 0 1 1;
    0 1 1 0;
    0 1 1 0;
    0 1 1 1;
    1 0 0 0;
    1 0 0 1;
    1 0 0 1;
    1 1 0 0;
    1 1 0 1;
    1 1 1 0;
    1 1 1 0;
    1 1 1 0;
    1 1 1 1;
    1 1 1 1;
    1 1 1 1];
xStr = char(reshape(1 : numel(x), size(x)));
for i = 1 : size(x, 1)
    xStr(i, :) = sprintf('%d', x(i, :));
end

NumRows = size(x, 1);
NumCols = size(x, 2);
prob = zeros(1, NumCols);
for i = 1 : length(prob)
    prob(i) = length(find(x(:, i))) / NumRows;
end

y = 0 : 2^(NumCols)-1;
yStr = dec2bin(y, NumCols);

probTrue = zeros(size(yStr, 1), 1);
for yrow = 1 : size(yStr, 1)
    for xrow = 1 : size(xStr, 1)
        if isequal(xStr(xrow, :), yStr(yrow, :))
            probTrue(yrow) = probTrue(yrow) + 1;
        end
    end
end
probTrue = probTrue / size(xStr, 1); % reproduce the middle column of Table I

probProduct = ones(size(yStr, 1), 1);
for row = 1 : size(yStr, 1)
    for col = 1 : size(yStr, 2)
        if yStr(row, col) == '0'
            probProduct(row) = probProduct(row) * (1 - prob(col));
        else
            probProduct(row) = probProduct(row) * prob(col);
        end
    end
end
display(probProduct) % reproduce the last column of Table I

ndx = find(probTrue ~= 0);
probProductError = sum(probTrue(ndx) .* log(probTrue(ndx) ./ probProduct(ndx)));
display(['Log-e error = ', num2str(probProductError)]); % 0.364, according to the end of Section IV
display(['Log-2 error = ', num2str(probProductError / log(2))]); 

% Calculate the mutual information between each pair of bits
% JointProb(xi, xj, vali, valj) = Prob(xi = vali-1, xj = valj-1)
Entropy = zeros(NumCols, NumCols);
JointProb = zeros(NumCols, NumCols, 2, 2);
for col1 = 1 : NumCols-1
    for col2 = col1+1 : NumCols
        for value1 = 0 : 1
            for value2 = 0 : 1
                for row = 1 : NumRows
                    if (x(row, col1) == value1) && (x(row, col2) == value2)
                        JointProb(col1, col2, value1+1, value2+1) = JointProb(col1, col2, value1+1, value2+1) + 1;
                    end
                end
            end
        end
        JointProb(col1, col2, :, :) = JointProb(col1, col2, :, :) / NumRows;
        for value1 = 0 : 1
            for value2 = 0 : 1
                JP = JointProb(col1, col2, value1+1, value2+1);
                if value1 == 0
                    Prob1 = 1 - prob(col1);
                else
                    Prob1 = prob(col1);
                end
                if value2 == 0
                    Prob2 = 1 - prob(col2);
                else
                    Prob2 = prob(col2);
                end
                Entropy(col1, col2) = Entropy(col1, col2) + JP * log(JP / Prob1 / Prob2);
            end
        end
    end
end
disp('Log-e entropy:')
disp(Entropy)
disp('Log-2 entropy:')
disp(Entropy / log(2))

% Calculate conditional probabilities
% CondProb(i, j, vi, vj) = Prob(xi = vi-1 | xj = vj-1)
% CondProb(i, i, vi, :)  = Prob(xi = vi-1)
CondProb = zeros(NumCols, NumCols, 2, 2);
for row = 1 : size(x, 1)
    for icol = 1 : size(x, 2)
        for jcol = 1 : size(x, 2)
            if icol == jcol
                CondProb(icol, icol, 1, :) = CondProb(icol, icol, 1, :) + (x(row, icol) == 0);
                CondProb(icol, icol, 2, :) = CondProb(icol, icol, 2, :) + (x(row, icol) == 1);
            else
                CondProb(icol, jcol, 1, 1) = CondProb(icol, jcol, 1, 1) + ((x(row, jcol) == 0) && (x(row, icol) == 0));
                CondProb(icol, jcol, 1, 2) = CondProb(icol, jcol, 1, 2) + ((x(row, jcol) == 1) && (x(row, icol) == 0));
                CondProb(icol, jcol, 2, 1) = CondProb(icol, jcol, 2, 1) + ((x(row, jcol) == 0) && (x(row, icol) == 1));
                CondProb(icol, jcol, 2, 2) = CondProb(icol, jcol, 2, 2) + ((x(row, jcol) == 1) && (x(row, icol) == 1));
            end
        end
    end
end
for icol = 1 : size(x, 2)
    for jcol = 1 : size(x, 2)
        if icol == jcol
            CondProb(icol, icol, :, :) = CondProb(icol, icol, :, :) / size(x, 1);
            disp(['Pr(x(', num2str(icol), ') = 0) = ', num2str(CondProb(icol, icol, 1, 1))]);
            disp(['Pr(x(', num2str(icol), ') = 1) = ', num2str(CondProb(icol, icol, 2, 1))]);
        else
            CondProb(icol, jcol, :, 1) = CondProb(icol, jcol, :, 1) ./ length(find(x(:, jcol) == 0));
            CondProb(icol, jcol, :, 2) = CondProb(icol, jcol, :, 2) ./ length(find(x(:, jcol) == 1));
            disp(['Pr(x(', num2str(icol), ') = 0 | x(', num2str(jcol), ') = 0) = ', num2str(CondProb(icol, jcol, 1, 1))]);
            disp(['Pr(x(', num2str(icol), ') = 1 | x(', num2str(jcol), ') = 0) = ', num2str(CondProb(icol, jcol, 2, 1))]);
            disp(['Pr(x(', num2str(icol), ') = 0 | x(', num2str(jcol), ') = 1) = ', num2str(CondProb(icol, jcol, 1, 2))]);
            disp(['Pr(x(', num2str(icol), ') = 1 | x(', num2str(jcol), ') = 1) = ', num2str(CondProb(icol, jcol, 2, 2))]);
        end
    end
end

probCol2 = ones(size(y))';
for i = 1 : length(y)
    bits = (yStr(i, :) == '1');
    probCol2(i) = CondProb(1, 1, bits(1)+1, 1) * CondProb(2, 1, bits(2)+1, bits(1)+1);
    probCol2(i) = probCol2(i) * CondProb(3, 2, bits(3)+1, bits(2)+1) * CondProb(4, 1, bits(4)+1, bits(1)+1);
end
disp('p(x1)*p(x2|x1)*p(x3|x2)*p(x4|x1):')
disp(probCol2) % reproduce the 2nd column of Table II
ndx = find(probTrue ~= 0);
probCol2Error = sum(probTrue(ndx) .* log(probTrue(ndx) ./ probCol2(ndx)));
display(['Log-e error = ', num2str(probCol2Error)]); % 0.094, according to the end of Section IV
display(['Log-2 error = ', num2str(probCol2Error / log(2))]);

probCol3 = ones(size(y))';
for i = 1 : length(y)
    bits = (yStr(i, :) == '1');
    probCol3(i) = CondProb(1, 1, bits(1)+1, 1) * CondProb(2, 1, bits(2)+1, bits(1)+1);
    probCol3(i) = probCol3(i) * CondProb(3, 2, bits(3)+1, bits(2)+1) * CondProb(4, 2, bits(4)+1, bits(2)+1);
end
disp('p(x1)*p(x2|x1)*p(x3|x2)*p(x4|x2):')
disp(probCol3) % reproduce the 3rd column of Table II
ndx = find(probTrue ~= 0);
probCol3Error = sum(probTrue(ndx) .* log(probTrue(ndx) ./ probCol3(ndx)));
display(['Log-e error = ', num2str(probCol3Error)]); % 0.094, according to the end of Section IV
display(['Log-2 error = ', num2str(probCol3Error / log(2))]);

probCol4 = ones(size(y))';
for i = 1 : length(y)
    bits = (yStr(i, :) == '1');
    probCol4(i) = CondProb(1, 1, bits(1)+1, 1) * CondProb(2, 1, bits(2)+1, bits(1)+1);
    probCol4(i) = probCol4(i) * CondProb(3, 2, bits(3)+1, bits(2)+1) * CondProb(4, 3, bits(4)+1, bits(3)+1);
end
disp('p(x1)*p(x2|x1)*p(x3|x2)*p(x4|x3):')
disp(probCol4) % reproduce the 4th column of Table II
ndx = find(probTrue ~= 0);
probCol4Error = sum(probTrue(ndx) .* log(probTrue(ndx) ./ probCol4(ndx)));
display(['Log-e error = ', num2str(probCol4Error)]); % 0.094, according to the end of Section IV
display(['Log-2 error = ', num2str(probCol4Error / log(2))]);

return