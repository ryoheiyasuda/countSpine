function [ point ] = findIntersection(highlightImageA,highlightImageB,xdim,ydim)
%FINDINTERSECTION Summary of this function goes here
%   Detailed explanation goes here
combination = highlightImageA + highlightImageB;
for i = 1:xdim
    for j = 1:ydim
        if combination(i,j) == 2 
            point = [i,j];
            return;
        end
    end
end
 

