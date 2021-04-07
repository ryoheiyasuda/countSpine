function a = cs_countSpines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters
smallSpine = 3; %eliminate spines smaller than this size.
longSpine = 30; %eliminate spines larger than this size
farSpine = 12; %pixel  Removes distant "spines" 
closeSpine = 4; %pixel Removes very small bumps
maxLength = 50; %pixel 30pixel = 
ZmaxLength = 4;
fw = 0;
fw2 = 2; %filter window for spines. (pixel)
resX = 0.5; %0.5; %resolution;
stuby_thin = 0.75; %micrometers;
dthresh = 0.25; %threshold of dendrite / dendritic intensity.
sthresh = 0.1; %threshold of spine intensity / dendritic intensity.
maxDistance = 1; %threshold for eliminating duplicates
debug = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Scale correction:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cs = get(gcf, 'UserData'); %Retrieves stored data
h1 = gcf; %gets the current figure handle
ha1 = gca; %gets the axes handle
filterWidnowSize = cs.param.filterWidnowSize; %Retrieves filter
sfactor = cs.param.sfactor; %Retrieves smoothing factor from cs structure
zfactor = cs.param.zfactor; 
FileName = cs.files.FileName; %Retrieves image file name
PathName = cs.files.PathName; % Retrieves image path name
mPerPixel = cs.param.mPerPixel; %Retrieves microns per pixel
mPerSlice = cs.param.mPerSlice; %Retrieves microns per z-slice
smallSpine = round(smallSpine*sfactor/resX); %eliminate spines smaller than this size.
farSpine = round(farSpine*sfactor/resX); %pixel
closeSpine = round(closeSpine*sfactor/resX);
longSpine = round(longSpine*sfactor/resX);

maxLength = round(maxLength*sfactor); %pixel
fw = round(fw*sfactor/resX); %filter window for adjacent curves.
fw2 = round(fw2*sfactor/resX);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Straighten the dendrite for simple calculation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Image2 = cs.ImageS;
%Image2(~BW) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    if strcmp(get(gco, 'String'), 'Count Spines') || strcmp(get(gco, 'String'), 'processing ...') %determines if "Count Spines" box is in either state
        hg = gco; %selects  "Count Spines" box
        %set(hg, 'Enable', 'Off');
        set(hg, 'String', 'processing ...'); %sets box to say "processing..."
        pause(0.1); %wait
    end
end
hs = get(gca, 'Children'); %gets all children of the axes handle(anything plotted in the figure on the axes)
for i=1:length(hs)
    tagstr = get(hs(i), 'Tag'); %gets the tag of each of the children of the axes handle
    if strcmp(tagstr, 'dendP')
        dendP = hs(i); %Unused
        dposy = get(hs(i), 'XData'); %Retrieves X Data from the previously drawn line
        dposx = get(hs(i), 'YData'); %Retrieves Y Data from the previously drawn line
    elseif strcmp(tagstr, 'NeuroImage')
        hI = hs(i); %Stores the handle for the image in hI
    elseif strcmp(tagstr, 'Spine') 
        delete(hs(i)); %Removes previously found spines(to repeat the process)
    elseif strcmp(tagstr, 'text')
        delete(hs(i));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rotation
spP2 = [dposx(1), dposy(1)]; %determines start/end point of line
spP1 = [dposx(end), dposy(end)]; %determines other start/end point of line
theta = -atan2(spP2(2)-spP1(2), spP2(1)-spP1(1)); %Determine the angle between the line drawn by points and the vertical 
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

    [cx, cy,p1] = improfile(max(Image2, [], 3), dposy + y1, dposx + x1, length(dposy)/resX); %Creates intensity profile crossing each point perpendicular to the start/end line
    im2 = [im2, p1(:)]; %Creates straigtened image
    xim = [xim, cx(:)];
    yim = [yim, cy(:)];
end

prf1 = ones(round(filterWidnowSize(1)/resX*2));
prf1 = prf1 / sum(prf1(:)); %Creates filter for image
im3= imfilter(im2, prf1, 'replicate'); %Creates final straigtened image
siz = size(im3); % determines dimensions of straightened image

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Background calculation based on the histogram.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:siz(2); 
    a1{i}=im2(:, i); 
    [n1,x1]=hist(a1{i}, max(a1{i})-min(a1{i}));  %Makes histograms of each vertical one-pixel wide segment
    [max1,pos]=max(n1); %Finds the mode of the histogram
    
    %b = mean(a1{i});
    try
        b(i) = x1(pos);  %%%%%%b is the baseline(mode of each vertical segment).
    catch
%       disp(['Error @ line 120: ', num2str(i)]); occurs frequently, does
        b(i) = 0;
    end
end
b = medfilt1(b, round(filterWidnowSize(1)/resX)); %adjusts baseline to account for image parameters
c = repmat(b(:)', [siz(1), 1]); %replicates 1-D baseline to the size of the image
im4 = im3 - c; %removes background from prefiltered, straigtened image

%im4 = im3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Counting spines using Fpeak program.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %round((filterWidnowSize(1)/resX + 1) / 8);
s = ceil(filterWidnowSize(1)/resX);
dtrim = s;
for i=1+fw:siz(2)-fw
    if mean(b(i-fw:i+fw)) < max(b)*dthresh %if this vertical segment is not part of the dendrite
        c = mean(im4(:, i-fw:i+fw), 2); %gets vertical segment
        c = imfilter(c, ones(1, s)/s, 'replicate');
        x = 1:length(c);
        try
            peak1{i}=fpeak(x,c,s,[dtrim,siz(1)-dtrim,max(b)*sthresh,1e10]); %outputs "peaks", places where the data has unusually high intensity
            peakx = peak1{i}(:,1); %temporarily stores peak locations
            peak1{i} = peak1{i}(peakx < siz(1)-dtrim & peakx > dtrim, :); %Removes peaks too close to the edge
        catch
%           disp(['error 148: ', num2str(i)]);
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

if debug    
    figure; imagesc(im3);
    hold on;
    for i=1:length(spineCell); 
        try
            plot(spineCell{i}(:,1), spineCell{i}(:,2), '-o', 'color', 'cyan', 'linewidth', 1);
        end
        try
            plot(spineCell2{i}(:,1), spineCell3{i}(:,2), '-o', 'color', 'white', 'linewidth', 2); 
        end
    end; %return;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Smoothing and drawing spines.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(spineCell3)
    spineCell3{i}(1:end-1,1) = imfilter(spineCell3{i}(1:end-1,1), prf, 'replicate');
    spineCell3{i}(1:end-1,2) = imfilter(spineCell3{i}(1:end-1,2), prf, 'replicate');
end

maxIm = max(Image2, [], 3);
set(hI, 'CData', maxIm);

%hold on; plot(dposy, dposx, '-', 'color', 'blue', 'linewidth', 2, 'Tag', 'Dendrite');
figure(h1);
axes(ha1);
mspine = 0;
sspine = 0;
tspine = 0;
for i=1:length(spineCell3)
        hold on;
        hspine(i) = plot(spineCell3{i}(:,1), spineCell3{i}(:,2), '-', 'color', 'white', 'linewidth', 2);
        spine_context = uicontextmenu;
        set(hspine(i), 'UIContextMenu', spine_context);
        uimenu(spine_context, 'Label', 'Stubby (red)', 'Callback', 'set(gco, ''color'', ''red''); cs_recalc');
        uimenu(spine_context, 'Label', 'Thin (green)', 'Callback', 'set(gco, ''color'', ''green''); cs_recalc');
        uimenu(spine_context, 'Label', 'Mushroom', 'Callback', 'set(gco, ''color'', ''white''); cs_recalc');
        uimenu(spine_context, 'Label', 'Delete', 'Callback', 'set(gco, ''color'', [0.5, 0.5, 0.5]); cs_recalc');

        dx = diff(spineCell3{i}(:,1));
        dy = diff(spineCell3{i}(:,2));
        r = sqrt(dx.^2 + dy.^2);
        length1 = sum(r);
        spineLength(i) = length1*mPerPixel;
        s_prof1 = spineCell3{i}(:,3);
        s_prof1 = s_prof1(1:round(length(s_prof1)*2/3));
        
        s_prof2  = improfile(maxIm, spineCell3{i}(1:end,1), spineCell3{i}(1:end,2), round(length1/resX));
        s_prof2 = imfilter(s_prof2, prf, 'replicate');
        %
        try
            peak1 = fpeak (1:length(s_prof2), s_prof2, fw2);
            maxval = peak1(1, 2);
            maxpos = peak1(1, 1); %%%% Pick the Far peak.
        catch
            [maxval, maxpos] = max(s_prof2);
        end

        spineInt2(i) = maxval;
        if maxpos > length(s_prof2)*0.6
            if spineLength(i) < stuby_thin
                spineType{i} = 'Stubby';
                set(hspine(i), 'color', 'red');
                sspine = sspine + 1;
            else
                spineType{i} = 'Thin';
                set(hspine(i), 'color', 'green');
                tspine = tspine + 1;
            end
        else
            spineType{i} = 'Mushroom';
            mspine = mspine + 1;
        end
        
        set(hspine(i), 'UserData', [spineInt1(i), spineInt2(i), spineLength(i)]);
        set(hspine(i), 'Tag', 'Spine');
end

% 

  
spineInt1 = spineInt1(:);
spineInt2 = spineInt2(:);
spineLength = spineLength(:);
spineNumber = length(spineInt1);

dx = diff(dposy);
dy = diff(dposx);
r = sqrt(dx.^2 + dy.^2);
dtrim = dtrim * resX; %Convert to pixel
dendLength = (sum(r)-2*dtrim)*mPerPixel;

a.spineInt1 = spineInt1;
a.spineInt2 = spineInt2;
a.spineType = spineType;
a.spineLength = spineLength;
a.averageIntensity = mean(spineInt2(:));
a.averageLength = mean(spineLength(:));
a.spineNumber = spineNumber;
a.mushroomSpine = mspine;
a.thinSpine = tspine;
a.stubbySpine = sspine;
a.dendLength = dendLength;
a.spineDensity = spineNumber / dendLength * 100;
a.spineCell = spineCell3;
 
% str{1} = ['FileName = ''', FileName, ''';'];
% str{2} = ['PathName = ''', PathName, ''';'];
% str{3} = ['mPerPixel = ', num2str(mPerPixel), ';'];
% str{4} = ['dtrim = ', num2str(mPerPixel), ';'];
cs.param.dtrim = dtrim;
cs.data = a;

set(h1, 'Toolbar', 'none');
set(h1, 'Toolbar', 'figure');
set(h1, 'UserData', cs);

uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                 'Position', [0.68, 0.0, 0.16, 0.04], 'String', 'Create Excel', 'Callback', ...
                 'spineExcel()', 'BackgroundColor', ...
                 [0.8,0.8,0.8]);           

cs_recalc;
try
    set(hg, 'Enable', 'On');
    set(hg, 'String', 'Count Spines');
end