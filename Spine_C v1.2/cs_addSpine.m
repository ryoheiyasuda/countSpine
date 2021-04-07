function cs_addSpine


cs = get(gcf, 'UserData');
FileName = cs.files.FileName;
PathName = cs.files.PathName;
mPerPixel = cs.param.mPerPixel;
mPerSlice = cs.param.mPerSlice;
dtrim = cs.param.dtrim;

[xdata, ydata, s_prof]  = improfile;
l = round(length(s_prof)*2/3);
spineInt1 = 0;
spineInt2 = max(s_prof(1:l));
dx = diff(xdata);
dy = diff(ydata);
r = sqrt(dx.^2 + dy.^2);
spineLength = sum(r)*mPerPixel;
    
hold on; 
hspine = plot(xdata, ydata, '--', 'color', 'white', 'linewidth', 2);
set(hspine, 'UserData', [spineInt1, spineInt2, spineLength]);
set(hspine, 'Tag', 'Spine');
spine_context = uicontextmenu;

set(hspine, 'UIContextMenu', spine_context);
uimenu(spine_context, 'Label', 'Stubby (red)', 'Callback', 'cs_set(gco,  ''red'')'); %original "set" and "delete" di not auto-recalc
uimenu(spine_context, 'Label', 'Thin (green)', 'Callback', 'cs_set(gco, ''green'')');
uimenu(spine_context, 'Label', 'Mushroom', 'Callback', 'cs_set(gco, ''white'')');
uimenu(spine_context, 'Label', 'Delete', 'Callback', 'cs_delete(gco)');

cs_recalc;