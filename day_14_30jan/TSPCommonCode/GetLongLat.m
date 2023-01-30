function [LongLat] = GetLongLat(Coord)
% Convert the DDD.MM format (where DDD is degrees and MM is minutes) of the data in the Coord array to degrees.
Deg = fix(Coord);
Min = Coord - Deg;
LongLat = pi * (Deg + 5 * Min / 3) / 180;
return