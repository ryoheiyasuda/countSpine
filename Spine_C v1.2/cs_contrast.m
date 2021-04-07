function cs_contrast(factor);

Clim1 = get(gca, 'Clim');
set(gca, 'Climmode', 'manual');
Clim1(2) = Clim1(2)/factor;
set(gca, 'Clim', Clim1);