function TestCountNormal

scaleFactor = 45.5*5; %45.5 um when zoom = 5;
filterWidnowSize = [3, 3]; %Ideally this should be smaller than PSF itself.
filtering = 'smooth'; %Smoothing often improves the quality of deconvolution.
%filtering = 'median';
%filtering = 'none';
part_img = 0; %If only a part of image is used
x_range = '1:256'; 
y_range = '1:256';
z_range = '5:32';

delpix = 3; %Pixels per segment.
maxLength = 30; %pixel 30pixel = 
ZmaxLength = 0; % Number of additional slides in stack to add to Max Projection
fw1 = 3; %Filterwindow for dendrite smooth.
backgroundAdjust = 200; %change for different microscopes (eliminates background over the mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reading file
[FileName,PathName] = uigetfile('*.tif','Select the tif-file'); %opens default user interface to open TIFF file
cd(PathName); 
fname = [PathName, FileName]; %sets filename
finfo1 = imfinfo([PathName, FileName]); %retrieves image
Zlen = length(finfo1); %determines stack size
Ylen = finfo1(1).Width; %Determines picture width
Xlen = finfo1(1).Height; %Determines picture height
try
    evalc(finfo1(1).ImageDescription); %attempts to read image parameters
    zoom = state.acq.zoomhundreds*100 + state.acq.zoomtens*10 + state.acq.zoomones; 
    mPerPixel = scaleFactor / zoom / state.acq.pixelsPerLine;  %micron per pixel
    mPerSlice = state.acq.zStepSize; %micron per slice
    rZ = mPerPixel / mPerSlice; %Ratio between lateral and axial pixels.
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recalibrate pixels.
    sfactor = zoom*state.acq.pixelsPerLine/5/512;
    delpix = round(delpix*sfactor); %Pixels per segment.
    maxLength = round(maxLength*sfactor); %pixel
    filterWidnowSize = round(filterWidnowSize*sfactor); %Ideally this should be smaller than PSF itself.
    zfactor = abs(1 / mPerSlice);
    %marzinz = round(marzinz * zfactor); %slices
    fw1 = round(fw1*sfactor); %pixelsFilterwindow for dendrite smooth.
    ZmaxLength = round(ZmaxLength * zfactor);

catch
    zoom = 3; %sets default parameters
    mPerPixel = 0.1;
    mPerSlice = 1;
    rZ = mPerPixel / mPerSlice;
    sfactor = 1;
    zfactor = 1;
end

%%%%%%%%%%%%%%%%%%%
cs.files.FileName = FileName; %Stores file name in cs structure
cs.files.PathName = PathName; %Stores path name in cs structure
cs.scaleFactor = scaleFactor; %Stores image scale in cs structure
cs.param.mPerPixel = mPerPixel; %Stores microns per pixel in cs structure
cs.param.mPerSlice = mPerSlice; %Stores microns per slice in cs structure
cs.param.sfactor = sfactor; %Stores smoothing factor in cs structure
cs.param.zfactor = zfactor; %Stores modification factor for z-slices in cs structure
cs.param.filterWidnowSize = filterWidnowSize; %Stores filter window size in cs structure
cs.param.delpix = delpix; %Stores pixels per segment in cs structure

% cs.param.smallSpine = smallSpine;
% cs.param.farSpine = farSpine;
% cs.param.closeSpine = closeSpine;
% cs.param.maxLength = maxLength;
% cs.param.filterWidnowSize = filterWidnowSize;
% cs.param.fw1 = fw1;
% cs.param.fw2 = fw2;
% cs.param.marzinz = marzinz;
% cs.param.ZmaxLength = ZmaxLength;

%%%%%%%%%%%%%%%%%%%
Image1 = zeros(Xlen, Ylen, Zlen, 'uint16'); %Creates matrix of zeroes 
Image2 = zeros(Xlen, Ylen, Zlen, 'uint16'); %Creates matrix of zeroes
for i=1:Zlen
    Image1(:,:,i) = imread(fname, i); %reads in each slice of TIFF file
    if strcmp (filtering, 'median') %reads in and filters each slice depending on given filter type
        Image2(:,:,i) = medfilt2(Image1(:,:,i), filterWidnowSize);
    elseif strcmp (filtering, 'smooth');
        Image2(:,:,i) = imfilter(Image1(:,:,i), ones(filterWidnowSize)/prod(filterWidnowSize), 'replicate');
    else
    end
end

cs.Image = Image1; %Stores original image in cs structure
cs.ImageF = Image2; %Stores filtered image in cs structure

bI = Image2(:, :,1); % first slice of filtered image
[n, x] = imhist(bI, 65535); %creates histogram of all pixels in first slice
[pos, background] = max(n); % finds the value of the mode(background intensity) of the image
background = uint16(background); %sets background to type uint16
%maxP = max(Image2(:)); %Unused
%threshold = (maxP - background)/12 + background; %Unused

if part_img %Determines if part of image is used
    evalc(['Image1 = Image1(', x_range, ',', y_range, ',', z_range, ')']);
    evalc(['Image2 = Image2(', x_range, ',', y_range, ',', z_range, ')']);
end

% 
background = background + backgroundAdjust; % increases background to include noise of greater intensity than the mode
Image2 = Image2 - background; %Subtracts background from image
Image2(Image2<0) = 0; %Sets all intensity values less than 0 in the image to 0
Image3 = Image2; %Creates copy of image

siz = size(Image2); %Determines dimensions of Image
Ylen = siz(1); Xlen = siz(2);  %Reverses width and height of image
if length(siz) > 2 %Sets Zlen equal to stack size
    Zlen = siz(3);
else
    Zlen = 1;
end
%%%%%%%%%%%%%%%%%
h1 = figure; %opens window
p1 = get(h1, 'position'); 
set(h1, 'position', [p1(1), p1(2) - p1(3) + p1(4), p1(3), p1(3)]); %sets place to put image
himage = imagesc(max(Image2,[],3)); %displays maxmimum projection
 
pause(0.1) %waits an instant
 waitforbuttonpress; % button down detected
 point1 = get(gca,'CurrentPoint');    
 waitforbuttonpress; % button down detected
 point2 = get(gca,'CurrentPoint');
 point1 = round(point1); %rounds point1
 point2 = round(point2); %rounds point2
 if point1(1) > point2(1) %Determines xpoint orientation
     maxx = point1(1) + maxLength; %finds zoomed window top point
     minx = point2(1) - maxLength; %finds zoomed window low point
 else
     minx = point1(1) - maxLength; %finds zoomed window low point
     maxx = point2(1) + maxLength; %finds zoomed window top point
 end
 if point1(3) > point2(3) %Determines ypoint orientation
     maxy = point1(3) + maxLength; %finds zoomed window right point
     miny = point2(3) - maxLength; %finds zoomed window left point
 else
     miny = point1(3) - maxLength; %finds zoomed window left point
     maxy = point2(3) + maxLength; %finds zoomed window right point
 end
% imZp1 = Image2(point1(1),point1(3),:);
% imZp2 = Image2(point2(1),point2(3),:);
% [q,p1Zmax] = max(imZp1(:));
% [q,p2Zmax] = max(imZp2(:));
% if p1Zmax > p2Zmax
%     maxz = p1Zmax + ZmaxLength;
%     minz = p2Zmax - ZmaxLength;
% else
%     maxz = p2Zmax + ZmaxLength;
%     minz = p1Zmax - ZmaxLength;
% end
if minx < 1; minx = 1; end; %sets low point to 1 if lower than 1
if miny < 1; miny = 1; end; %sets left point to 1 if lower than 1
% if minz < 1; minz = 1; end;
if maxx > Xlen; maxx = Xlen; end %sets high point to zoomed image ceiling if higher than zoomed image ceiling
if maxy > Ylen; maxy = Ylen; end %sets high point to right cutoff if higher than right cutoff
Image2Store = Image2; %Stores Image2 temporairily 
Image3Store = Image3; %Stores Image3 temporairily  
Image2 = Image2(miny:maxy, minx:maxx, 1:end); %crops Image2
Image3 = Image3(miny:maxy, minx:maxx, 1:end); %crops Image3
maxProj = max(Image3, [], 3); %creates maximum projection of Image3
figure(h1); %brings up window
hI = imagesc(maxProj); %Displays cropped window 
zoomxSize = maxx - minx; %finds zoomed image height 
zoomySize = maxy - miny; %finds zoomed image width
if ((point1(1) < point2(1) && point1(3) < point2(3)) || (point1(1) > point2(1) && point1(3) > point2(3))) %Determines proper point orientation
origpoint1 = [maxLength,maxLength]; %Adjusts first clicked point for new zoomed image
origpoint2 = [zoomxSize - maxLength, zoomySize - maxLength];  %Adjusts second clicked point for new zoomed image
else
    origpoint1 = [maxLength, zoomySize - maxLength]; %Adjusts first clicked point for new zoomed image
    origpoint2 = [zoomxSize - maxLength, maxLength]; %Adjusts second clicked point for new zoomed image
end
[xi, yi] = findMid(maxProj,origpoint1,origpoint2); %Determines pathway
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Strip image for faster calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Image2 = Image2Store;
Image3 = Image3Store; %Returns Image3 to original Image
%imagesc(max(Image3, [], 3));
for i=1:length(xi) 
    posx = round(yi(i)); %rounds each y value
    posy = round(xi(i)); %rounds each x value
    imZ = Image2(posx, posy, :); %takes one-dimensional slice along Z-dimension
    [v, posz] = max(imZ(:)); % determines index of maximum intensity along z-dimension
    sp = [posx, posy, posz];  %stores point temporarily
    dend_pos{i} = sp; %Creates cell for each point
    dposx(i) = sp(1); %stores x value of each point
    dposy(i) = sp(2); %stores y value of each point
    dposz(i) = sp(3); %stores z value of each point
end
% 

%error = 0;
dend_pos_end = sp; %Determines final point (value left in "sp")
%%%
%count = 1; %Unused
%sp = round(dend_pos{1}); %Unused
%sp2 = round(dend_pos_end); %Unused

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Repeats image crop process for ALL points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dposz = zFilter(dposz); %Remove outliers (3D top-down intersections)
xs = min([dposx-maxLength]);
xe = max([dposx+maxLength]);
ys = min([dposy-maxLength]);
ye = max([dposy+maxLength]);
zs = min([dposz-ZmaxLength]);
ze = max([dposz+ZmaxLength]);
if xs < 1; xs = 1; end;
if ys < 1; ys = 1; end;
if zs < 1; zs = 1; end;
if xe > Xlen; xe = Xlen; end
if ye > Ylen; ye = Ylen; end
if ze > Zlen; ze = Zlen; end

%zs = 1; ze = Zlen;

Image2 = Image2(1:end, 1:end, zs:ze);
%Image3 = Image3(xs:xe, ys:ye, zs:ze);

% siz = size(Image2);
% Xlen = siz(1); Ylen = siz(2); Zlen = siz(3);
% [X, Y, Z] = meshgrid(1:Ylen, 1:Xlen, 1:Zlen);
% xr = X-sp(2); yr = Y-sp(1); zr = Z-sp(3);
% R = round(sqrt(xr.^2 + yr.^2));

for i=1:length(dend_pos)
    dend_pos{i} = dend_pos{i} - [xs, ys, zs] + [1, 1, 1];
end
dposx = dposx - xs + 1;
dposy = dposy - ys + 1;
%dposz = dposz - zs + 1;

cs.ImageS = Image2; %Saves Image2 for counting

%%%%%Creates User Interface%%%%%

%figure(h1); 
set(hI, 'Tag', 'NeuroImage');
%plot(dposy, dposx, '-', 'color', 'white');
radius = 1;
try
[dposx,dposy] = filterRepeats(dposx,dposy); %Note: minor bug, does not take inputs for identical X values, very rare, fix
[yy, xx]=cs_spline(dposy, dposx);
catch
    disp(dposx);
    disp(dposy);
end
cs.figure.dendP = line(yy, xx, 'color', 'white', 'Tag', 'dendP', 'linewidth', 2);
step = round(length(yy)/15);
yy1 = yy(1:step:end);
xx1 = xx(1:step:end);
% if mod(length(xx), step)
%     yy1 = [yy1, yy(end)];
%     xx1 = [xx1, xx(end)];
% end

%Creates dragable polyline

for i=1:length(xx1)
   roiPos = [yy1(i)-radius, xx1(i)-radius, radius*2, radius*2];
   cs.figure.dendH = rectangle('Position', roiPos, 'Tag', 'dendH', 'EdgeColor', 'black', 'FaceColor', 'White', ...
       'ButtonDownFcn', 'cs_lineDrag', 'Curvature', [1,1]);
end



set(h1, 'UserData', cs);

uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                'Position', [0.0, 0.0, 0.25, 0.08], 'String', 'Count Spines', 'Callback', ...
                'cs_countSpines', 'BackgroundColor', [0.8,0.8,0.8]); 
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                'Position', [0.84, 0.0, 0.16, 0.04], 'String', 'Recalc', 'Callback', ...
                'cs_recalc', 'BackgroundColor', [0.8,0.8,0.8]); 
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                'Position', [0.84, 0.04, 0.16, 0.04], 'String', 'Add spine', 'Callback', ...
                'cs_addSpine', 'BackgroundColor', [0.8,0.8,0.8]);            
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                 'Position', [0.92, 0.12, 0.08, 0.04], 'String', 'Length', 'Callback', ...
                 'cs_recalc(''Length'')', 'BackgroundColor', [0.8,0.8,0.8]);
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                 'Position', [0.92, 0.16, 0.08, 0.04], 'String', 'Intensity', 'Callback', ...
                 'cs_recalc(''Intensity'')', 'BackgroundColor', ...
                 [0.8,0.8,0.8]);
             
uicontrol ('Style', 'text', 'Unit', 'normalized', ...
                 'Position', [0.92, 0.5, 0.08, 0.04], 'String', 'Contrast', 'Callback', ...
                 '', 'BackgroundColor', ...
                 [0.8,0.8,0.8]);
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                 'Position', [0.92, 0.46, 0.08, 0.04], 'String', 'Up', 'Callback', ...
                 'cs_contrast(2)', 'BackgroundColor', ...
                 [0.8,0.8,0.8]);
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                 'Position', [0.92, 0.42, 0.08, 0.04], 'String', 'Down', 'Callback', ...
                 'cs_contrast(1/2)', 'BackgroundColor', ...
                 [0.8,0.8,0.8]);   
uicontrol ('Style', 'pushbutton', 'Unit', 'normalized', ...
                'Position', [0.68, 0.04, 0.16, 0.04], 'String', 'Length Calculator', 'Callback', ...
                'lengthcalc', 'BackgroundColor', [0.8,0.8,0.8]);           
             
return;
