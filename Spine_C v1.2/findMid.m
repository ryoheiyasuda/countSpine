function [ yi, xi ] = findMid( maxImage,origpoint1,origpoint2 )
%FINDPATH Summary of this function goes here
%   Detailed explanation goes here
 maxImage1 = imresize(maxImage,0.5);
 point1 = round(origpoint1*0.5);
 point2 = round(origpoint2*0.5);
[xdim,ydim] = size(maxImage1);
realImageThresh = 1200;
maxImage1 = maxImage1 - realImageThresh;
highlightImageA = zeros(xdim, ydim);
highlightImageA = treatSides(highlightImageA,xdim,ydim);
highlightImageB = highlightImageA;
highlightImageA(point1(2),point1(1)) = 1;
highlightImageB(point2(2),point2(1)) = 1;
maxImage1(point1(2),point1(1)) = 0;
maxImage1(point2(2),point2(1)) = 0;
[xclass, yclass, largediff] = classify(point1, point2);
debugcount = 1;
maxImage1 = treatSidesZero(maxImage1,xdim,ydim);
maxImage2 = maxImage1; 
positionImageA = cell(xdim,ydim);
positionImageA(:) = {NaN};
positionImageB = positionImageA;
while noIntersection(highlightImageA,highlightImageB)
[highlightImageA, maxImage1, positionImageA] = highlightMaxSurroundings(maxImage1, highlightImageA, positionImageA, xclass, yclass,  1, largediff, xdim, ydim);
[highlightImageB, maxImage2, positionImageB] = highlightMaxSurroundings(maxImage2, highlightImageB, positionImageB, xclass, yclass,  0, largediff, xdim, ydim);
highlightImageA(highlightImageA > 1) = 1;
highlightImageB(highlightImageB > 1) = 1;
debugcount = debugcount + 1;
if debugcount == 200
    return;
end
end
midpoint = findIntersection(highlightImageA,highlightImageB,xdim,ydim);
path = findPath(midpoint, positionImageA, positionImageB);
for i = 1 : length(path)
        xi(i) = round(path{i}(1)*2);
        yi(i) = round(path{i}(2)*2);
end
end


   

    