% figure properties for computational bottleneck paper

% format = 'png'; % for panels tricky for EPS (e.g. Pcolor plots)
% color = 'rgb';
% dpi = 600;

fontsize = 9;
fontname = 'Arial';
M = 3; % marker size for univariate scatter plots
sym = 'o';  % markers for scatters and strip plots

Units = 'centimeters';

% line widths
widths.plot = 0.75;
widths.error = 0.5;
widths.axis = 0.5;


% panel sizes
figsize.square = [4 4];
figsize.smallmatrices = [3 3];
figsize.colourbar_height = 0.05;

% colours for matrix plots
colourmaps.gamma_error = brewermap(5,'Blues');
colourmaps.trials = brewermap(10,'BuPu');  % colormap for all gamma except gamma=1
colourmaps.stability = brewermap(10,'OrRd');  % colormap for all gamma except gamma=1

colourmaps.no_decay = brewermap(10,'Greys');


% colours for lines
colourmaps.gamma_lines = brewermap(6,'*Set1');

% exportpath
exportpath = '../Panels/';
