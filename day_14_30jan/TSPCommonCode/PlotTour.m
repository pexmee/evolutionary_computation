function PlotTour(Coord, Tour)
% Plot a TSP tour. 
% INPUTS: Coord = n x 2 array of lat/long coordinates of n cities
%         Tour = n-element ordered array of cities to visit
figure, hold on
for i = 1 : length(Tour)-1
    plot(Coord(Tour(i),1), Coord(Tour(i),2), 'r*');
    plot([Coord(Tour(i),1) Coord(Tour(i+1),1)], [Coord(Tour(i),2), Coord(Tour(i+1),2)], 'b-');
end
xlabel('Longitude')
ylabel('Latitude')
box on
return