function [] = recurFind(point1, point2, maxImage, i1, midi, i2)
%RECURDRAW Summary of this function goes here
%   Detailed explanation goes here
global points
    if i2 - midi == 1
        points{midi} = findMid(points{i1}, points{i2}, maxImage);
    return
    else
        points{midi} = findMid(points{i1}, points{i2}, maxImage);
        recurFind(point1, points{midi}, maxImage, i1, (i1 + midi) / 2, midi);
        recurFind(points{midi}, point2, maxImage, midi, (midi + i2) / 2, i2);
    end
end

