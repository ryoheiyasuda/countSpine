function [spineInt1, spineInt2, spineLength, spineNumber, spineCell3] = countSpines;
LR = 0.4; %Larteral resolution in micrometers;
AR = 1; %Axial resolution; %The resolution is changed by an addative algorithm.
scaleFactor = 83*3; %83 um when zoom = 3;
delpix = 2; %Pixels per segment.
smallSpine = 3; %eliminate spines smaller than this size.
farSpine = 15; %pixel

%maxImageSize = 256; %%If the image exceeds this size, the program calculates each 256x256 segment separately.
filtering = 'smooth'; %Smoothing often improves the quality of deconvolution.
%filtering = 'median';
%filtering = 'none';
filterWidnowSize = [3, 3]; %Ideally this should be smaller than PSF itself.
fw1 = 4; %Filterwindow for dendrite smooth.

part_img = 1; %If only a part of image is used
x_range = '1:256'; 
y_range = '1:256';
z_range = '5:32';

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading file
[FileName,PathName] = uigetfile('*.tif','Select the tif-file');
finfo1 = imfinfo([PathName, '/', FileName]);
Zlen = length(finfo1);
Ylen = finfo1(1).Width;
Xlen = finfo1(1).Height;
evalc(finfo1(1).ImageDescription);
zoom = state.acq.zoomhundreds*100 + state.acq.zoomtens*10 + state.acq.zoomones;
mPerPixel = scaleFactor / zoom / state.acq.pixelsPerLine;  %micron per pixel
mPerSlice = state.acq.zStepSize; %micron per slice
rZ = mPerPixel / mPerSlice; %Ratio between lateral and axial pixels.

%%%%%%%%%%%%%%%%%%%
Image1 = zeros(Xlen, Ylen, Zlen, 'uint16');
Image2 = zeros(Xlen, Ylen, Zlen, 'uint16');
for i=1:Zlen
    Image1(:,:,i) = imread([PathName, '/', FileName], i); %-background;
    if strcmp (filtering, 'median')
        Image2(:,:,i) = medfilt2(Image1(:,:,i), filterWidnowSize);
    elseif strcmp (filtering, 'smooth');
        Image2(:,:,i) = imfilter(Image1(:,:,i), ones(filterWidnowSize)/prod(filterWidnowSize), 'replicate');
    else
    end
end

bI = Image2(:, :,1);
[n, x] = imhist(bI, 65535);
[pos, background] = max(n);
background = uint16(background);
maxP = max(Image2(:));
threshold = (maxP - background)/12 + background;

if part_img
    evalc(['Image1 = Image1(', x_range, ',', y_range, ',', z_range, ')']);
    evalc(['Image2 = Image2(', x_range, ',', y_range, ',', z_range, ')']);
end


siz = size(Image2);
Ylen = siz(1); Xlen = siz(2); Zlen = siz(3);

% 
Image2 = Image2 - background;
Image2(Image2<0) = 0;
Image3 = Image2;

siz = size(Image2);
Ylen = siz(1); Xlen = siz(2); Zlen = siz(3);

h1 = figure; imagesc(max(Image2,[],3)); % creates initial color image
pause(0.1)
waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
point1 = point1(1,1:2);    
point1 = round(point1);
posx = point1(2);
posy = point1(1);
imZ = Image2(posx, posy, :);
[v, posz] = max(imZ(:));
sp = [posx, posy, posz]; %Gives (x,y,z) value of starting point
count = 1;
dend_pos{count} = sp;
dposx(count) = sp(1);
dposy(count) = sp(2);
dposz(count) = sp(3);
[X, Y, Z] = meshgrid(1:Xlen, 1:Ylen, 1:Zlen); 
xr = X-sp(2); yr = Y-sp(1); zr = Z-sp(3);
R = round(sqrt(xr.^2 + yr.^2));
Image2(R <= delpix & abs(zr) <= 1) = 0;
error = 0;

waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    % button down detected
point1 = point1(1,1:2);    
point1 = round(point1);
posx = point1(2);
posy = point1(1);
imZ = Image2(posx, posy, :);
[v, posz] = max(imZ(:));
sp = [posx, posy, posz];
count = 2;
dend_pos{count} = sp;
dposx(count) = sp(1);
dposy(count) = sp(2);
dposz(count) = sp(3);
xr = X-sp(2); yr = Y-sp(1); zr = Z-sp(3); % confirm
R = round(sqrt(xr.^2 + yr.^2));        %confirm
Image2(R <= delpix & abs(zr) <= 1) = 0;
error = 0;
%sp = dend_pos{1};

while count < 100 & ~error
    imA = Image2; %(xr, yr, zr);
    xr = X-sp(2); yr = Y-sp(1); zr = Z-sp(3);
    R = round(sqrt(xr.^2 + yr.^2));
    imA( (R > delpix+1 | R < delpix) | (abs(zr) > 1)) = 0;
    spP2 = dend_pos{count-1};
    spP1 = dend_pos{count};
    theta = -atan2(spP2(2)-spP1(2), spP2(1)-spP1(1));
    x1 = cos(theta)*xr + sin(theta)*yr;
    y1 = -sin(theta)*xr + cos(theta)*yr;
    imA(y1 > 0) = 0;
    [val, pos] = max(imA(:));
    if val ~= 0
        posz = floor ((pos + Xlen*Ylen - 0.5) / Xlen / Ylen);
        posxy = round(pos - (posz-1) * Xlen*Ylen);
        posy = floor ((posxy + Xlen - 0.5) / Xlen);
        posx = round(posxy - (posy-1)*Xlen);
        count = count + 1;
        dend_pos{count} = [posx, posy, posz];
        sp = dend_pos{count};
        dposx(count) = sp(1);
        dposy(count) = sp(2);
        dposz(count) = sp(3);
        Image2(R <= delpix & abs(zr) <= 1) = 0;
        if sp(1)<= delpix*2 | sp(1) >= Xlen - delpix*2 | ...
                sp(2) <= delpix*2 | sp(2) >= Ylen - delpix*2 | ...
                sp(3) <= 2 | sp(3) >= Zlen - 2
            error = 1;
        end
        
    else
        error =1 ;
    end
end
dposx = dposx(3:end);
dposy = dposy(3:end);
dposz = dposz(3:end);

prf = ones(fw1);
prf = prf / sum(prf(:));
dposx = imfilter(dposx, prf, 'replicate');
dposy = imfilter(dposy, prf, 'replicate');
dposz = imfilter(dposz, prf, 'replicate');

dend_pos = {};
count = count - 2;
for i=1:count
    dend_pos{i} = [dposx(i), dposy(i), dposz(i)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Finish dendrite drawing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Image2 = Image3;
%Rotation
spP2 = dend_pos{1};
spP1 = dend_pos{end};
theta = -atan2(spP2(2)-spP1(2), spP2(1)-spP1(1));

im2 = [];
xim = [];
yim = [];
resX = 0.5;
for i=-25:resX:25;
    xr = 0;
    yr = i;
    x1 = cos(theta)*xr + sin(theta)*yr;
    y1 = -sin(theta)*xr + cos(theta)*yr;

    [cx, cy,p1] = improfile(max(Image2, [], 3), dposy + y1, dposx + x1, length(dposy)*delpix/resX);
    im2 = [im2, p1(:)];
    xim = [xim, cx(:)];
    yim = [yim, cy(:)];
end

im3= imfilter(im2, ones(filterWidnowSize*2)/prod(filterWidnowSize), 'replicate');
siz = size(im3);
im4 = [];

for i=1:siz(2); 
    a1{i}=im2(:, i); 
    [n1,x1]=hist(a1{i}, max(a1{i})-min(a1{i})); 
    [max1,pos]=max(n1); 
    b(i) = x1(pos); 
    c = a1{i} - b(i);
    im4=[im4, c];

end

fw = round((filterWidnowSize(1)/resX + 1) / 4);
s = filterWidnowSize(1)/resX;
for i=1+fw:siz(2)-fw
    if mean(b(i-fw:i+fw)) < max(b)* 0.15
        c = mean(im4(:, i-fw:i+fw), 2);
        x = 1:length(c);
        try
            peak1{i}=fpeak(x,c,s,[fw,siz(1)-fw,max(b)/25,1e10]);
        catch
            disp(i);
            figure; plot(x, c);
            peak1{i} = [];
        end
    else
        peak1{i} = [];
    end
end

s = filterWidnowSize(1)/resX;
hsegment = s;
vsegment = s*2;
spineCell = segment_spines(peak1, vsegment, hsegment);
bh = max(b)/2;
bwidth = (1 + length(find(b > bh)))/2;
bwidth = round(bwidth);
cent = round((length(b) + 1)/2);

%figure; imagesc(im3); 

spine = 0;
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
    pixval = spineCell{i}(:,3); 
    aa = find(pixval >= max(pixval(:))/2);
    spineStart = aa(1);
    if length(pixval(spineStart:end)) <= smallSpine
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
        
        for j = 1:length(x1)
            spineCell3{spine}(j, 1) = xim(x1(j), y1(j));
            spineCell3{spine}(j, 2) = yim(x1(j), y1(j));
            spineCell3{spine}(j, 3) = v1(j);
            %spineCell3{spine}(j, 3) = spineCell2{spine}(j,3);
        end
    end
end

for i=1:length(spineCell3)
    spineCell3{i}(1:end-1,1) = imfilter(spineCell3{i}(1:end-1,1), prf, 'replicate');
    spineCell3{i}(1:end-1,2) = imfilter(spineCell3{i}(1:end-1,2), prf, 'replicate');
end

maxIm = max(Image2, [], 3);
figure (h1); 
hold on; plot(dposy, dposx, '-', 'color', 'blue', 'linewidth', 2);
for i=1:length(spineCell3)
        hold on;
        plot(spineCell3{i}(:,1), spineCell3{i}(:,2), '-', 'color', 'white');
        s_prof  = improfile(maxIm, spineCell3{i}(1:end-1,1), spineCell3{i}(1:end-1,2));
       
        spineInt1(i) = max(spineCell3{i}(:,3));
        spineInt2(i) = max(s_prof);
        dx = diff(spineCell3{i}(:,1));
        dy = diff(spineCell3{i}(:,2));
        r = sqrt(dx.^2 + dy.^2);
        spineLength(i) = sum(r)*mPerPixel;
end
spineInt1 = spineInt1(:);
spineInt2 = spineInt2(:);
spineLength = spineLength(:);
spineNumber = length(spineInt1);