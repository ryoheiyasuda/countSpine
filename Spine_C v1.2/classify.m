function [ xclass, yclass, largediff ] = classify( point1, point2 )
%CLASSIFY Summary of this function goes here
%   Detailed explanation goes here
yresidual = point2(1) - point1(1);
xresidual = point2(2) - point1(2);

if xresidual > 0 
    xclass = 1;
else
    xclass = 0;
end

if yresidual > 0
    yclass = 1;
else
    yclass = 0;
end

if xresidual > yresidual
    largediff = 1;
else
    largediff = 0;
end

