function [xi, yi] = findAndDrawPoints( maxImage,origpoint1,origpoint2 )
%FINDPOINTS Summary of this function goes here
%   Detailed explanation goes here
 maxImage1 = imresize(maxImage,0.5);
 origpoint1 = round(origpoint1*0.5);
 origpoint2 = round(origpoint2*0.5);
[yi, xi] = findMid(origpoint1, origpoint2, maxImage1);
  %  [xi,yi] = adjustPoints(xiPre,yiPre,maxImage);
  %  dendLength =  (83.7/512) * otherPolylineLength(points1); %change value for different microscopes (pixel/um ratio)
  %  str1 = sprintf('Dendrite Length :  %3.1f um', dendLength);
  %  display =  text (5, xdim-5, str1, 'color', 'white', 'VerticalAlignment', 'bottom');
  %  set(display, 'Interpreter', 'none');
end

