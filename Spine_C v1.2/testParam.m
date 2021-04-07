function [ spineNum ] = testParam( testFile, testdp, testsp, testmaxp )
%TESTPARAM Summary of this function goes here
%   Detailed explanation goes here

%Parameters
load(testFile);
dthresh = testdp; %threshold of dendrite / dendritic intensity.
sthresh = testsp; %threshold of spine intensity / dendritic intensity.
maxDistance = testmaxp; %threshold for eliminating duplicates

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%try
%    if strcmp(get(gco, 'String'), 'Count Spines') | strcmp(get(gco, 'String'), 'processing ...') 
%        hg = gco;
        %set(hg, 'Enable', 'Off');
%        set(hg, 'String', 'processing ...');
%        pause(0.1);
%    end
% end
%hs = get(gca, 'Children');
%for i=1:length(hs)
%    tagstr = get(hs(i), 'Tag');
%    if strcmp(tagstr, 'dendP')
%        dendP = hs(i);
%        dposy = get(hs(i), 'XData');
%        dposx = get(hs(i), 'YData');
%    elseif strcmp(tagstr, 'NeuroImage')
%        hI = hs(i);
%    elseif strcmp(tagstr, 'Spine')
%        delete(hs(i));
%    elseif strcmp(tagstr, 'text')
%        delete(hs(i));
%    end
%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotation
spP2 = [dposx(1), dposy(1)];
spP1 = [dposx(end), dposy(end)];
theta = -atan2(spP2(2)-spP1(2), spP2(1)-spP1(1));
% startz = floor(min(dposz)-marzinz);
% endz = ceil(max(dposz)+marzinz);
% if startz < 1
%     startz = 1;
% end
% if endz > Zlen
%     endz = Zlen;
% end
% Image2 = Image2(:,:,startz:endz);

im2 = [];
xim = [];
yim = [];

for i=-maxLength:resX:maxLength;
    xr = 0;
    yr = i;
    x1 = cos(theta)*xr + sin(theta)*yr;
    y1 = -sin(theta)*xr + cos(theta)*yr;

    [cx, cy,p1] = improfile(max(Image2, [], 3), dposy + y1, dposx + x1, length(dposy)/resX);
    im2 = [im2, p1(:)];
    xim = [xim, cx(:)];
    yim = [yim, cy(:)];
end

prf1 = ones(round(filterWidnowSize(1)/resX*2));
prf1 = prf1 / sum(prf1(:));
im3= imfilter(im2, prf1, 'replicate');
siz = size(im3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Background calculation based on the histogram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:siz(2); 
    a1{i}=im2(:, i); 
    [n1,x1]=hist(a1{i}, max(a1{i})-min(a1{i})); 
    [max1,pos]=max(n1);
    
    %b = mean(a1{i});
    try
        b(i) = x1(pos);  %%%%%%b is the baseline.
    catch
        disp(['Error @ line 120: ', num2str(i)]);
        b(i) = 0;
    end
end
b = medfilt1(b, round(filterWidnowSize(1)/resX));
c = repmat(b(:)', [siz(1), 1]);
im4 = im3 - c;

%im4 = im3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Counting spine. using Fpeak program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %round((filterWidnowSize(1)/resX + 1) / 8);
s = ceil(filterWidnowSize(1)/resX);
dtrim = s;
for i=1+fw:siz(2)-fw
    if mean(b(i-fw:i+fw)) < max(b)*dthresh
        c = mean(im4(:, i-fw:i+fw), 2);
        c = imfilter(c, ones(1, s)/s, 'replicate');
        x = 1:length(c);
        try
            peak1{i}=fpeak(x,c,s,[dtrim,siz(1)-dtrim,max(b)*sthresh,1e10]);
            peakx = peak1{i}(:,1);
            peak1{i} = peak1{i}(peakx < siz(1)-dtrim & peakx > dtrim, :);
        catch
            disp(['error 148: ', num2str(i)]);
%             figure; plot(x, c);
%             peak1{i} = [];
        end
    else
        peak1{i} = [];
    end
end

%s = filterWidnowSize(1)/resX;
xsegment = s; %ceil(s/2);  %x_difference.
if xsegment < 2;
    xsegment = 2;
end
ysegment = s*2;
spineCell = segment_spines(peak1, xsegment, ysegment, closeSpine);


%return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Removing "wrong" spines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bh = max(b)/2;
bwidth = (1 + length(find(b > bh)))/2;
bwidth = round(bwidth);
bw2 = (1 + length(find (b > max(b)*sthresh)))/2;
bw2 = round(bw2);
cent = round((length(b) + 1)/2);

spine = 0;
prf = ones(fw2);
prf = prf/sum(prf(:));

for i=1:length(spineCell);
    siz = size(spineCell{i});
    badSpine = 0;
    if siz(1) <= smallSpine %Remove very small spines.
        badSpine = 1;
    end
    lastPixelY = spineCell{i}(end,2);
    if spineCell{i}(end, 1) < cent
        lastPixelX = cent - bwidth;
        if spineCell{i}(end, 1) < cent - bwidth - farSpine
            badSpine = 1;
        end
    else
        lastPixelX = cent + bwidth;
        if spineCell{i}(end, 1) > cent + bwidth + farSpine
            badSpine = 1;
        end
    end
    if spineCell{i}(1, 1) <=3 | spineCell{i}(1, 1) >= length(b)-3
        badSpine = 1;
    end
    x1 = spineCell{i}(:,1);
    y1 = spineCell{i}(:,2);
    pixval = spineCell{i}(:,3);
    try
        peak1 = fpeak([1:length(pixval)], pixval, fw2);
        maxval = peak1(1,2);
        maxvals = peak1(:,2);
        if length(maxvals) >= 2
            maxmaxval = max(maxvals);
            maxvals = maxvals(maxvals > maxmaxval/5);
            maxval = maxvals(1);
        end
    catch
        maxval = max(pixval);
        %disp(i);
    end
    aa = find(pixval >= maxval/2);
    spineStart = aa(1);
    if length(pixval(spineStart:end)) <= smallSpine
        badSpine = 1;
    end
    if x1(spineStart) < cent - bwidth - longSpine | x1(spineStart) > cent + bwidth + longSpine
        badSpine = 1;
    end
    
    if x1(spineStart) > cent - bwidth - closeSpine & x1(spineStart) <   cent
        badSpine = 1;
    elseif x1(spineStart) < cent + bwidth + closeSpine & x1(spineStart) > cent
        badSpine = 1;
    end
    lastPixelV = 0;
    if ~badSpine
        spine = spine+1;
        spineCell2{spine} = spineCell{i};     
        spineCell2{spine} = [spineCell2{spine}; lastPixelX, lastPixelY, 0];

        x1 = spineCell2{spine}(spineStart:end,2);
        y1 = spineCell2{spine}(spineStart:end,1);
        v1 = spineCell2{spine}(spineStart:end,3);

        %hold on;
        %plot(y1, x1, 'color', 'white');
        spineInt1(spine) = maxval;
        for j = 1:length(x1)
            spineCell3{spine}(j, 1) = xim(x1(j), y1(j));
            spineCell3{spine}(j, 2) = yim(x1(j), y1(j));
            spineCell3{spine}(j, 3) = v1(j);
            %spineCell3{spine}(j, 3) = spineCell2{spine}(j,3);
        end
    end
end

spineCell3 = distanceFilter(spineCell3, maxDistance); %removes duplicate spines based on proximity

spineNum = length(spineCell3);

