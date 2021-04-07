function [betterLines] = distanceFilter(allLines, maxDistance)
count = 1;
ignorenum = 0 ;
ignore = [];
for g = 1 : (length(allLines) - 1)
    if any(ignore  == g)
        continue
    end
    for h = (g + 1) : length(allLines)
         if any(ignore == h)
             continue
         end 
        currentDistance = 100000000;
        for i = 1 : length(allLines{g})
            for j = 1 : length(allLines{h})
                 newDistance = sqrt(((allLines{g}(i,1) - allLines{h}(j,1))^2) + ((allLines{g}(i,2) - allLines{h}(j,2))^2));
                    if  currentDistance > newDistance
                        currentDistance = newDistance;
          
                    end
            end
        end
        if (currentDistance < maxDistance)
            if sqrt(((allLines{g}(end,1) - allLines{g}(1,1))^2) + ((allLines{g}(end,2) - allLines{g}(1,2))^2)) > ...
                sqrt(((allLines{h}(end,1) - allLines{h}(1,1))^2) + ((allLines{h}(end,2) - allLines{h}(1,2))^2))
                ignorenum = ignorenum + 1;
                ignore(ignorenum) = h; 
            else
                ignorenum = ignorenum + 1;
                ignore(ignorenum) = g; 
            end
        end
    end
end
betterLines = allLines;
betterLines(ignore) = [];
        
           
        