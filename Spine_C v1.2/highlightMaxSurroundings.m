function [newHighlightImage, newMaxImage, newPositionImage] = highlightMaxSurroundings(maxImage,highlightImage, positionImage, xclass, yclass,  num, ld, xdim, ydim)

%%%%%%Determine "for" format%%%%%%%

if xclass + num == 1
    xcl = xdim:-1:1;
else
    xcl = 1:xdim;
end
if yclass + num == 1 
    ycl = ydim:-1:1;
else
    ycl = 1:ydim;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ld
for i = xcl
    for j = ycl
        currentPoint = highlightImage(i,j);
        if currentPoint == 1
            [maxX,maxY] = findMaxPix(maxImage, i,j);
            if ~isnan(maxX) && ~isnan(highlightImage(maxX,maxY)) && (~(maxX ==  i) || ~(maxY == j)) 
                highlightImage(maxX,maxY) = 2;
                positionImage(maxX,maxY) = {[i,j]};
                maxImage(maxX,maxY) = 0;
            elseif ~isnan(maxX) && ~isnan(highlightImage(maxX,maxY))
                highlightImage(maxX,maxY) = 2;
                maxImage(maxX,maxY) = 0;
            else
                highlightImage(maxX,maxY) = NaN;
            end
        end
    end
end
else
    for i = ycl
        for j = xcl
            currentPoint = highlightImage(j,i);
            if currentPoint == 1
                [maxX,maxY] = findMaxPix(maxImage, j,i);
                if ~isnan(maxX) && ~isnan(highlightImage(maxX,maxY)) && (~(maxY ==  i) || ~(maxX == j)) 
                    highlightImage(maxX,maxY) = 2;
                    positionImage(maxX,maxY) = {[j,i]};
                    maxImage(maxX,maxY) = 0;
                elseif ~isnan(maxX) && ~isnan(highlightImage(maxX,maxY))
                    highlightImage(maxX,maxY) = 2;
                    maxImage(maxX,maxY) = 0;
                else
                    highlightImage(maxX,maxY) = NaN;
                end
            end
        end
    end
end



newPositionImage = positionImage;
newHighlightImage = highlightImage;
newMaxImage = maxImage;



