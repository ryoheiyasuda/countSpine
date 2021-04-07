function [ newImage ] = treatSides( image, xdim, ydim )
%TREATSIDES Summary of this function goes here
%   Detailed explanation goes here
for i = 1:xdim
    for j = 1:ydim
        if i == 1 || j == 1 || i == xdim || j == ydim
            image(i,j) = NaN;
        end
    end
end
newImage = image;


