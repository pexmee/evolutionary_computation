function PlotBestTour(Filename, Coord, DistanceArray)
% Plot the best tour from a .TOUR file, with format defined in http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95.
% INPUTS: Filename = name of file (no extension)
%         Coord = n x 2 array of lat/long coordinates of n cities
%         DistanceArray = n x n array of distances between cities
fid = fopen([Filename, '.opt.tour']);
FileData = textscan(fid, '%s');
fclose(fid);
FileData = FileData{1, 1};
for i = 1 : length(FileData)
    if isequal(FileData{i}, 'TOUR_SECTION'), break, end
end
NumCities = str2double(FileData{i-1});
Tour = zeros(1, NumCities+1);
j = i + 1;
for i = 1 : NumCities
    Tour(i) = str2double(FileData{j});
    j = j + 1;
end
Tour(end) = Tour(1);
Distance = CalcDistance(Tour, DistanceArray);
disp(['TSPLIB95 Minimum distance = ', num2str(Distance)]);
PlotTour(Coord, Tour), title('TSPLIB95 Best Tour')
return