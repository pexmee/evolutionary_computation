function [DistanceArray] = CreateDistanceArray(Coord, EdgeWeightType)
% Given the n x 2 array Coord containing the latitude/longitude coordinates of n cities,
% calculate the n x n array of distances between each pair of cities.
% Input EdgeWeightType indicates what distance calculation to use
dim = size(Coord, 1);
DistanceArray = zeros(dim, dim);
RRR = 6378.388;
for i = 1 : dim
    for j = 1 : dim
        if isequal(EdgeWeightType, 'EUC_2D')
            DistanceArray(i, j) = round( norm(Coord(i,:) - Coord(j,:)) ); % two-dimensional Euclidean distance
        elseif isequal(EdgeWeightType, 'GEO')
            q1 = cos( Coord(i,2) - Coord(j,2) );
            q2 = cos( Coord(i,1) - Coord(j,1) );
            q3 = cos( Coord(i,1) + Coord(j,1) );
            DistanceArray(i, j) = floor( RRR * acos( 0.5*((1+q1)*q2 - (1-q1)*q3) ) + 1 ); % Geographical distance
        end
    end
end
return