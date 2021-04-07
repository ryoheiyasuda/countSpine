function [ filteredZList ] = zFilter( zList )
%ZFILTER Summary of this function goes here
%   Detailed explanation goes here
 startingPoint = zList(1);
 filteredZList = zList(1);
 for i = 2:length(zList)
     if (abs(startingPoint - zList(i)) < 7)
         startingPoint = zList(i);
         filteredZList = [filteredZList,zList(i)];         
     end
 end
end

