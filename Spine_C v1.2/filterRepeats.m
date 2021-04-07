function [ newxx, newyy ] = filterRepeats( xx,yy )
%FILTERREPEATS Summary of this function goes here
%   Detailed explanation goes here
for i = 1:length(xx)
xy{i} = [xx(i),yy(i)];
end
for i = 1 : (length(xx) - 1)
    for  j = i + 1 : length(xx)
        if xy{i} == xy{j}
            [newxx,newyy] = filterRepeats([xx(1:j-1),xx(j+1:end)],[yy(1:j-1),yy(j+1:end)]);
            return;
        end
    end
end
newxx = xx;
newyy = yy;

end

