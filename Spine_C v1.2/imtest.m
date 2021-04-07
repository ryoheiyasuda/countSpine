function [ newIm ] = imtest( Image )
%IMTEST Summary of this function goes here
%   Detailed explanation goes here
for i = 1 :  size(Image,3)
    newIm(:,:,i)  = imadjust(Image(:,:,i));
end

    

end

