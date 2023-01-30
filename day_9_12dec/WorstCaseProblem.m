function [InitFunction, CostFunction, FeasibleFunction] = WorstCaseProblem
InitFunction = @Init;
CostFunction = @Cost;
FeasibleFunction = [];
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population, OPTIONS] = Init(OPTIONS)
% Initialize population
Population = struct('chrom', cell([1 OPTIONS.popsize]), 'cost', cell([1 OPTIONS.popsize]));
for popindex = 1 : OPTIONS.popsize
    Population(popindex).chrom = randi([0, 1], [1, OPTIONS.numVar]);
    if OPTIONS.Gray
        Population(popindex).chrom = BinaryToGray(Population(popindex).chrom);
    end
end
OPTIONS.OrderDependent = false;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Population] = Cost(Population, OPTIONS)
% Compute the cost of each member in Population
for popindex = 1 : length(Population)
    gene = BitsToGene(Population(popindex).chrom, OPTIONS);
    if mod(gene, 2) == 0
        Population(popindex).cost = 1; % cost for even numbers
    else
        Population(popindex).cost = 2; % cost for odd numbers
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gene] = BitsToGene(Bits, OPTIONS)
% Convert the bits in the Bits array to a real-valued gene, assuming binary coding
if OPTIONS.Gray
    Bits = GrayToBinary(Bits);
end
Gene = sum(Bits .* (2.^(0 : length(Bits)-1)));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Binary] = GrayToBinary(Gray)
% Convert the gray-coded bit string in Gray to a binary-coded bit string in Binary.
% The first bit in each array is the LSB.
Binary = zeros(size(Gray));
Binary(end) = Gray(end);
for i = length(Gray)-1 : -1 : 1
    Binary(i) = mod(Binary(i+1) + Gray(i), 2);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Gray] = BinaryToGray(Binary)
% Convert the binary-coded bit string in Binary to a gray-coded bit string in Gray.
% The first bit in each array is the LSB.
Gray = Binary;
for i = length(Gray)-1 : -1 : 1
    Gray(i) = xor(Binary(i+1), Binary(i));
end
return