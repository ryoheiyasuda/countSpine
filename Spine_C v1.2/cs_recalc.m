function cs_recalc (str)

if ~nargin
    str = '';
end

cs = get(gcf, 'UserData');
FileName = cs.files.FileName;
PathName = cs.files.PathName;
mPerPixel = cs.param.mPerPixel;
mPerSlice = cs.param.mPerSlice;
dtrim = cs.param.dtrim;
c = get(gca, 'Children');

j = 0;
mspine = 0;
sspine = 0;
tspine = 0;
for i=1:length(c)
    tagS = get(c(i), 'Tag');
    typeS = get(c(i), 'Type');
    if strcmp(tagS, 'Spine') & strcmp(typeS, 'line')
        uData = get(c(i), 'UserData');
        XData = get(c(i), 'XData');
        YData = get(c(i), 'YData');
        isspine = 1;
        if sum(get(c(i), 'color') == [1,1,1]) == 3
            spineT = 'Mushroom';
            mspine = mspine + 1;
        elseif sum(get(c(i), 'color') == [0,1,0]) == 3 %Blue
            spineT = 'Thin';
            tspine = tspine + 1;
        elseif sum(get(c(i), 'color') == [1,0,0]) == 3 %Red
            spineT = 'Stubby';
            sspine = sspine + 1;
        elseif sum(get(c(i), 'color') == [0.5, 0.5, 0.5]) == 3 %gray
            isspine = 0;
        end
        if isspine
            j = j+1;
            spineInt1(j) = uData(1);
            spineInt2(j) = uData(2);
            spineLength(j) = uData(3);
            spineType{j} = spineT;
            if XData(1) < XData(end)
                x1 = XData(1) - 4;
            else
                x1 = XData(1);
            end
            if YData(1) < YData(end)
                y1 = YData(1)-2;
            else
                y1 = YData(1)+2;
            end
            text(x1, y1, num2str(j), 'color', 'white');
        end
    elseif strcmp(tagS, 'dendP')
        dposx = get(c(i), 'XData');
        dposy = get(c(i), 'YData');
    elseif strcmp(get(c(i), 'Type'), 'text')
        delete(c(i));
    elseif strcmp(get(c(i), 'Type'), 'NeuroImage')
    end
end

siz = size(cs.ImageS);
Ylen = siz(1);

spineInt1 = spineInt1(:);
spineInt2 = spineInt2(:);
spineLength = spineLength(:);
spineNumber = length(spineInt1);

dx = diff(dposy);
dy = diff(dposx);
r = sqrt(dx.^2 + dy.^2);
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
a.spineCell = cs.data.spineCell;
h1 = gcf;
cs.data = a;
set(h1, 'UserData', cs);
str1 = sprintf('Density (#spines/100um): %3.1f (%d spines / %3.1f um)', a.spineDensity, a.spineNumber, a.dendLength);
str2 = ['Mushroom: ', num2str(a.mushroomSpine), '   Stubby (red): ', num2str(a.stubbySpine), '   Thin (green): ', num2str(a.thinSpine)];
str3 = sprintf('Mean spine intensity %3.1f', a.averageIntensity);
str4 = sprintf('Mean spine length %3.2f um', a.averageLength);
str5 = sprintf('%s\n%s\n%s\n%s\n%s\n%s', str1, str2, str3, str4, PathName, FileName);
blind = 'Blind Mode';
ht = text (5, Ylen-5, blind, 'color', 'white', 'VerticalAlignment', 'bottom');
set(ht, 'Interpreter', 'none');
%disp (str5);

if strcmp (str, 'Intensity')
    disp('Intensity');
    for i=1:length(spineInt2)
        str6 = sprintf('%8.2f', spineInt2(i));
        disp(str6);
    end
elseif strcmp (str, 'Length')
    disp('Length (um)');
    for i=1:length(spineInt2)
        str6 = sprintf('%8.2f', spineLength(i));
        disp(str6);
    end
end