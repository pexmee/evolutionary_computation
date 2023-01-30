function [Distance] = CalcDistance(Tour, DistanceArray)
% Calculate the distance of a TSP tour.
% INPUTS: Tour = n-element ordered vector of city indices to visit on the tour
%         DistanceArray = n x n matrix of distances between cities
Distance = 0;
for i = 1 : length(Tour)-1
    Distance = Distance + DistanceArray(Tour(i), Tour(i+1));
end
return