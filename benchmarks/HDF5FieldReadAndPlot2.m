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
HDF5_path = 'D:\TEMP\OCCIGEN_CUPRAD\compares4\';
HDF5_filename = "results_GfP.h5";


kr = 1;

Nx_plot = 2;
Ny_plot = 5;

rhorange = [0 0.3e-3];

rhoplot_range = 1:40:200;
kz_plot = 5;

print_A4 = true;

resol='-r300';
colors={'-g','-k','--r'};
set(0,'defaulttextInterpreter','latex')

%% computed params
HDF5_filepath = strcat(HDF5_path,HDF5_filename);



%% data



Fields_rzt = h5read(HDF5_filepath,"/IRprop/Fields_rzt");
output_field = h5read(HDF5_filepath,"/outputs/output_field");
output_plasma = h5read(HDF5_filepath,"/outputs/output_plasma");
StartFieldsR = h5read(HDF5_filepath,"/pre-processed/startfield_r");
tgrid = h5read(HDF5_filepath,"/IRprop/tgrid");
rgrid = h5read(HDF5_filepath,"/IRprop/rgrid");
zgrid = h5read(HDF5_filepath,"/IRprop/zgrid");

Nz = length(zgrid); Nt = length(tgrid); Nr = length(rgrid);


%% tests
firstplasma = squeeze(output_plasma(1,:,:));
firstplasma(1,end)

figure
pcolor(firstplasma)
shading interp
colorbar


%% plot first_intens
fig.sf(1).method = @plot;
fig.sf(1).arg{1} = tgrid;
fig.sf(1).arg{2} = squeeze(abs(output_field(1,kr,:)));
I_tmp = fig.sf(1).arg{2};

fig.sf(2) = fig.sf(1);
t_1e = 1e-15*h5read(HDF5_filepath,"/inputs/laser_pulse_duration_in_1_e");
fig.sf(2).arg{2} = max(I_tmp)*exp(-(tgrid/t_1e).^2);

fig.title = 'First intensity';
plot_preset_figure(fig,'default');


return

%% plot z_evol
fig.sf(1).method = @plot;
fig.sf(1).arg{1} = tgrid;
fig.sf(1).arg{2} = squeeze(output_field(1,kr,:));

for k1 = 2:Nz
    fig.sf(k1) = fig.sf(1);
    fig.sf(k1).arg{2} = squeeze(output_field(k1,kr,:));
end


fig.title = 'z-evol';
plot_preset_figure(fig,'default');


%% plot r_evol
clear fig;
fig.sf(1).method = @plot;
fig.sf(1).arg{1} = tgrid;
fig.sf(1).arg{2} = squeeze(output_field(kz_plot,rhoplot_range(1),:));
fig.sf(1).legendlabel = strcat('$\rho = ', num2str(1e6*rgrid(1)), '~\mathrm{\mu m}$');

for k1 = rhoplot_range(2:end)
    fig.sf(k1) = fig.sf(1);
    fig.sf(k1).arg{2} = squeeze(output_field(kz_plot,k1,:));
    fig.sf(k1).legendlabel = strcat('$\rho = ', num2str(1e6*rgrid(k1)), '~\mathrm{\mu m}$');
end


fig.title = 'r-evol';
plot_preset_figure(fig,'default');


%% plot startfield
clear fig;
fig.sf(1).method = @pcolor; fig.sf(1).shading = 'interp';
fig.sf(1).arg{1} = tgrid;
fig.sf(1).arg{2} = rgrid;
fig.sf(1).arg{3} = StartFieldsR';

fig.title = 'startfield real';
plot_preset_figure(fig,'default');

%% plot z-full-field
clear fig;
fig.sf(1).method = @pcolor; fig.sf(1).shading = 'interp';
fig.sf(1).arg{1} = tgrid;
fig.sf(1).arg{2} = rgrid;
fig.sf(1).arg{3} = squeeze(output_field(kz_plot,:,:));

fig.title = 'chosen field real';
plot_preset_figure(fig,'default');

if print_A4

%% print t-fields

k_plot = 1;
k2 = 1;
for k1 = 1:Nz
    arr_fig.fig(k2).sf(1).method = @plot; 
    arr_fig.fig(k2).sf(1).arg{1} = tgrid; arr_fig.fig(k2).sf(1).arg{2} = squeeze(output_field(1,kr,:));
    if ((mod(k1,Nx_plot*Ny_plot)==0) || (k1 == Nz) )
        arr_fig.filenamepng = strcat('Frzt_',num2str(k_plot),'.png');
        arr_fig.resolutionpng = '-r450';
        Print_Array_Figs_A4(arr_fig,2,5);
        k_plot = k_plot + 1;
        k2 = 0;
        clear arr_fig;
    end
    k2 = k2+1;
end

k_plot = 1;
k2 = 1;
for k1 = 1:Nz
    arr_fig.fig(k2).sf(1).method = @plot; 
    arr_fig.fig(k2).sf(1).arg{1} = tgrid; arr_fig.fig(k2).sf(1).arg{2} = squeeze(Fields_rzt(1,kr,:));
    if ((mod(k1,Nx_plot*Ny_plot)==0) || (k1 == Nz) )
        arr_fig.filenamepng = strcat('Frzt_',num2str(k_plot),'_check.png');
        arr_fig.resolutionpng = '-r450';
        Print_Array_Figs_A4(arr_fig,2,5);
        k_plot = k_plot + 1;
        k2 = 0;
        clear arr_fig;
    end
    k2 = k2+1;
end

%% print r-fields
k_plot = 1;
k2 = 1;
for k1 = 1:Nz
    arr_fig.fig(k2).sf(1).method = @plot; arr_fig.fig(k2).xlim = rhorange;
    arr_fig.fig(k2).sf(1).arg{1} = rgrid; arr_fig.fig(k2).sf(1).arg{2} = squeeze(output_field(k1,:,ceil(Nt/2)));
    if ((mod(k1,Nx_plot*Ny_plot)==0) || (k1 == Nz) )
        arr_fig.filenamepng = strcat('rfield_',num2str(k_plot),'.png');
        arr_fig.resolutionpng = '-r450';
        Print_Array_Figs_A4(arr_fig,2,5);
        k_plot = k_plot + 1;
        k2 = 0;
        clear arr_fig;
    end
    k2 = k2+1;
end

%% print plasma
k_plot = 1;
k2 = 1;
for k1 = 1:Nz
    arr_fig.fig(k2).sf(1).method = @pcolor; arr_fig.fig(k2).sf(1).shading = 'interp'; arr_fig.fig(k2).sf(1).colorbar = 'eastoutside';
    %arr_fig.fig(k2).ylim = [-rhorange(2), rhorange(2)];
%     arr_fig.fig(k2).sf(1).arg{1} = tgrid; arr_fig.fig(k2).sf(1).arg{2} = rgrid; arr_fig.fig(k2).sf(1).arg{3} = squeeze(output_plasma(k1,:,:));
    arr_fig.fig(k2).sf(1).arg{1} = 1e6*squeeze(output_plasma(k1,:,:)); 
    if ((mod(k1,Nx_plot*Ny_plot)==0) || (k1 == Nz) )
        arr_fig.filenamepng = strcat('plasma_',num2str(k_plot),'.png');
        arr_fig.resolutionpng = '-r450';
        Print_Array_Figs_A4(arr_fig,2,5);
        k_plot = k_plot + 1;
        k2 = 0;
        clear arr_fig;
    end
    k2 = k2+1;
end


%% print inst_intens
k_plot = 1;
k2 = 1;
for k1 = 1:Nz
    arr_fig.fig(k2).sf(1).method = @pcolor; arr_fig.fig(k2).sf(1).shading = 'interp'; arr_fig.fig(k2).sf(1).colorbar = 'eastoutside';
    %arr_fig.fig(k2).ylim = [-rhorange(2), rhorange(2)];
%     arr_fig.fig(k2).sf(1).arg{1} = tgrid; arr_fig.fig(k2).sf(1).arg{2} = rgrid; arr_fig.fig(k2).sf(1).arg{3} = abs(squeeze(output_field(k1,:,:))).^2;
    arr_fig.fig(k2).sf(1).arg{1} = (abs(squeeze(output_field(k1,:,:))).^2)';
    if ((mod(k1,Nx_plot*Ny_plot)==0) || (k1 == Nz) )
        arr_fig.filenamepng = strcat('inst_intens_',num2str(k_plot),'.png');
        arr_fig.resolutionpng = '-r450';
        Print_Array_Figs_A4(arr_fig,2,5);
        k_plot = k_plot + 1;
        k2 = 0;
        clear arr_fig;
    end
    k2 = k2+1;
end


end

return
%%
kz = 1;
kr = 1;

Erz1i = squeeze(StartFieldsR(kr,:) );

Erz1ir = squeeze(StartFieldsR(:,512) );
% Erz2i = squeeze(StartFieldsR(kr+5,:) );

Erz1 = squeeze(StartFieldsR(1,:,1));

Erz1r = squeeze(StartFieldsR(:,512,1));

figure
plot(tgrid,Erz1i);  hold on
% plot(tgrid,Erz2i);  hold on

figure
plot(rgrid,Erz1ir);  hold on

figure
plot(tgrid,Erz1);  hold on

figure
plot(rgrid,Erz1r);  hold on


return

Erz1 = squeeze( Fields_rzt(kz,kr,:) );
Erz2 = squeeze( Fields_rzt(kz+1,kr,:) );

figure
plot(tgrid,Erz1);  hold on
% plot(tgrid,Erz2)

return

% plot(tgrid,Erz1,'-*');  hold on
% plot(tgrid,Erz2,'-*')

% plot(Erz1,'-*');  hold on
% plot(Erz2,'-*')

rgrid(kr)
zgrid(kz)


kz = 1;
kt = 1024;

Ezt1 = squeeze( Fields_rzt(kz,:,kt) );
Ezt2 = squeeze( Fields_rzt(kz+1,:,kt) );
Ezt3= squeeze( Fields_rzt(kz+2,:,kt) );

figure
plot(rgrid,Ezt1); hold on
% plot(rgrid,Ezt2); hold on
% plot(rgrid,Ezt3); hold on

% plot(rgrid,Ezt1,'-*'); hold on
% plot(rgrid,Ezt2,'-*'); hold on
% plot(rgrid,Ezt3,'-*'); hold on

% plot(Ezt1,'-*'); hold on
% plot(Ezt2,'-*'); hold on
% plot(Ezt3,'-*'); hold on

tgrid(kt)
zgrid(kz)


%%%%%%%%%%%%%



% Splot = squeeze(S(:,:,k_zplot));
% pcolor(Hgrid,rgrid,abs(Splot))
% colorbar
% shading interp
% colormap jet

% caxis([-38 -20])
% caxis([0 1e-10])
% caxis([0 1e-5])
% ylim([0 1000*0.02])
% caxis([0 3e-11])

% %     xlim([0 5]);
% %     ylim([-80 80])
%     xlabel('Harmonic order [-]')
%     ylabel('$\varphi$ [mrad] (taken in $D = 1~\mathrm{m}$)')
%     title('$|\mathcal{E}|$ [arb.u]')
% %     print('SpectraOnScreen',resol,'-dpng')
 
    
    
    
% figure
% Szplot = squeeze(S(:,k_oplot,:));
% % pcolor(zgrid,rgrid,abs(Szplot))
% pcolor(zgrid,rgrid,log(abs(Szplot)))
% colorbar
% shading interp
% colormap jet

return
    
    xlim([14 22]);
%     ylim([-80 80])
    print('SpectraOnScreenZoom',resol,'-dpng')

colorbar