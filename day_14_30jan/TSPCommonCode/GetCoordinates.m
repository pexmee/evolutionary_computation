function [Coord, EdgeWeightType] = GetCoordinates(Filename)
% Get the TSP coordinates from a .TSP file, with format defined in http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95.
% This code is not completely general - it fails for some of the .TSP files from TSPLIB95 because there are small differences
% in format between some of the files.
% OUTPUTS: Coord = n x 2 array of lat/long coordinates of cities
%          EdgeWeightType = string defining how to calculate distances between cities
fid = fopen([Filename, '.tsp']);
FileData = textscan(fid, '%s');
fclose(fid);
FileData = FileData{1, 1};
for i = 1 : length(FileData)
    if isequal(FileData{i}, 'DIMENSION:'), break, end
end
dim = str2double(FileData{i+1});
for k = i+2 : length(FileData)
    if isequal(FileData{k}, 'EDGE_WEIGHT_TYPE:'), break, end
end
EdgeWeightType = FileData{k+1};
Coord = zeros(dim, 2);
for j = k+2 : length(FileData)
    if isequal(FileData{j}, 'NODE_COORD_SECTION'), break, end
end
for i = 1 : dim
    Coord(i, 1) = str2double(FileData{j+2});
    Coord(i, 2) = str2double(FileData{j+3});
    j = j + 3;
end
return