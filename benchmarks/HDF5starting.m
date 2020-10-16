clear all; close all;
% clearvars

%% head
addpath('D:\Dropbox\PhD\Analyses\shared_functions') 
currfold = pwd;

global Ip_HeV hbar alpha_fine c_light elcharge elmass rBohr mu0 eps0;
Ip_HeV = 27.21138602; hbar = 1.054571800e-34; alpha_fine=1/137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass = 9.10938356e-31;
mu0=4*pi*1e-7; eps0=1/(mu0*c_light^2);
rBohr = 4*pi*eps0*hbar^2/(elmass*elcharge^2);

global EFS INTENSITYau
EFS = 5.1422e11; % 1*a.u. = EFS*V/m
INTENSITYau = (1/(8*pi*alpha_fine))*hbar^3/(elmass^2*rBohr^6);


delete('*.png');

%% parameters

% path = 'D:\data\CUPRAD';
HDF5_path = 'D:\TEMP\OCCIGEN_CUPRAD\compares2\';
HDF5_filename = "results_short.h5";


kr = 2;

Nx_plot = 2;
Ny_plot = 5;

rhorange = [0 0.3e-3];

rhoplot_range = 1:40:200;
kz_plot = 2;

print_A4 = true;

resol='-r300';
colors={'-g','-k','--r'};
set(0,'defaulttextInterpreter','latex')

%% computed params
HDF5_filepath = strcat(HDF5_path,HDF5_filename);



%% data



StartFieldsR = h5read(HDF5_filepath,"/pre-processed/startfield_r");





%% plot startfield
clear fig;
fig.sf(1).method = @pcolor; fig.sf(1).shading = 'interp';
% fig.sf(1).arg{1} = tgrid;
% fig.sf(1).arg{2} = rgrid;
fig.sf(1).arg{1} = StartFieldsR;

fig.title = 'startfield real';
plot_preset_figure(fig,'default');


%% plot cut
clear fig;
fig.sf(1).method = @plot;
% fig.sf(1).arg{1} = tgrid;
% fig.sf(1).arg{2} = rgrid;
fig.sf(1).arg{1} = StartFieldsR(kr,:);

fig.title = 'startfield real cut';
plot_preset_figure(fig,'default');


