function [fit] = genFit(presult, ideal)
fit = 0;
for 1 : length(presult)
    fit = fit + ((presult - ideal)^2);
end

