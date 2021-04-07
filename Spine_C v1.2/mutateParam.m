function [ newD, newS, newMax ] = mutateParam(dthresh, sthresh, maxDistance)
r1 = (round(rand(1)*2) - 1);
r2 = (round(rand(1)*2) - 1);
r3 = (round(rand(1)*2) - 1);
newD = dthresh;
newS = sthresh;
newMax = maxDistance;
while ~(r1 == 0)
newD = newD + (0.01 * r1);
r1 = (round(rand(1)*2) - 1);
end
while ~(r2 == 0)
newS = newS + (0.01 * r2);
r2 = (round(rand(1)*2) - 1);
end
while ~(r3 == 0)
newMax = newMax + (0.1 * r3);
r3 = (round(rand(1)*2) - 1);
end

