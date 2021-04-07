function [plength] = polylineLength(polyline)
plength = 0;
for i = 1 : (length(polyline) - 1)
    plength = plength + sqrt(((polyline(i,1) - polyline((i + 1),1))^2) + ((polyline(i,2) - polyline((i+1),2))^2));
end

