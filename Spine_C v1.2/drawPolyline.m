function [ points ] = drawPolyline( points )
%DRAWLINE Summary of this function goes here
%   Detailed explanation goes here

for i = 1 : length(points)
    xcoords(i) = points{i}(1);
    ycoords(i) = points{i}(2);
end
hold on; plot(xcoords, ycoords, '-', 'color', 'blue', 'linewidth', 2);

end
