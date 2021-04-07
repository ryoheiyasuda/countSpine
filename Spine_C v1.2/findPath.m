function [ path ] = findPath(mid, posA, posB )
%FINDPATH Summary of this function goes here
%   Detailed explanation goes here
pathBack = {};
pathForward = {};
currentPoint = mid;
while ~isnan(posA{currentPoint(1),currentPoint(2)})
    currentPoint = posA{currentPoint(1),currentPoint(2)};
    pathBack = [currentPoint,pathBack];   
end
currentPoint = mid;
while ~isnan(posB{currentPoint(1),currentPoint(2)})
    currentPoint = posB{currentPoint(1),currentPoint(2)};
    pathForward = [pathForward,currentPoint];
end
path = [pathBack, mid, pathForward];
end

