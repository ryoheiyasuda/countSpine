function [  ] = lengthcalc(  )
%LENGTHCALC Summary of this function goes here
%   Detailed explanation goes here
pixToumConversionFactor = 0.1432;
[x,y] = ginput(2);
 distance = sqrt(((x(1) - x(2))^2) + ((y(1) - y(2))^2))*pixToumConversionFactor;
 comment = ['The measured length is ' num2str(distance) ' um.'];
disp(comment);
end

