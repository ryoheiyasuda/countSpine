function [ maxX, maxY ] = findMaxPix(maxImage, x1,y1)
count = 1;
maxInt = 0;
xadd = 0; %instantiated because of NaN inequality properties
yadd = 0; %instantiated because of NaN inequality properties
for i = -1:1:1
    for j = -1:1:1
        if maxImage(x1 + i, y1 + j) > maxInt
            maxInt = maxImage(x1 + i, y1 + j);
            xadd = i;
            yadd = j;
        end
    end
end
        if (maxInt < 1) 
            maxX = x1;
            maxY = y1;
            return;
        else 
            maxX = x1 + xadd;
            maxY = y1 + yadd;
            return;
        end
            
            


