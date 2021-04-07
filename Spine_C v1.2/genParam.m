function [ finald, finals, finalMax ] = genParam( numGens )
%GENPARAM Summary of this function goes here
%   Detailed explanation goes here
dp = 0.18;
sp = 0.13;
maxp = 0.7;
childrenNum = 10; %number of children per generation
ideal = [47,29,18]; % ideal solution set
for h = 1 : numGens
    for i = 1 : childrenNum  
        for j = 1:3
            resultList(j) = testParam(['pic', num2str(j), 'data.mat'],dp,sp,maxp);
        end
        mutFit = 0;
        for k = 1:3
            mutFit = mutFit + ((resultList(k) - ideal(k))^2);
        end
        dpList(i) = dp;
        spList(i) = sp;
        maxpList(i) = maxp;
        mutFitList(i) = mutFit;
        [dp,sp,maxp] = mutateParam(dp,sp,maxp);  
    end 
    [bestFit, bestFiti]  = min(mutFitList);
    dp = dpList(bestFiti);
    sp = spList(bestFiti);
    maxp = maxpList(bestFiti);
end
finald = dp;
finals = sp;
finalMax = maxp;

