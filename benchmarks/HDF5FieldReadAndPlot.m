clear all; close all;

global Ip_HeV hbar alpha_fine c_light elcharge elmass;
Ip_HeV = 27.21138602; hbar = 1.054571800e-34; alpha_fine=1/137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass = 9.10938356e-31;

global EFS
EFS = 5.1422e11; % 1*a.u. = EFS*V/m


% path = 'D:\data\CUPRAD';
path = 'D:\TEMP\ECLIPSE_CUPRAD';

resol='-r300';
colors={'-g','-k','--r'};
set(0,'defaulttextInterpreter','latex')

%% plot

currfold = pwd;
cd(path)
% AllFields = h5read("results.h5","/IRprop/Fields_rzt");
% tgrid = h5read("results.h5","/IRprop/tgrid");
% rgrid = h5read("results.h5","/IRprop/rgrid");
% zgrid = h5read("results.h5","/IRprop/zgrid");

% AllFields = h5read("results_reference.h5","/IRprop/Fields_rzt");
% tgrid = h5read("results_reference.h5","/IRprop/tgrid");
% rgrid = h5read("results_reference.h5","/IRprop/rgrid");
% zgrid = h5read("results_reference.h5","/IRprop/zgrid");

AllFields = h5read("results4.h5","/IRprop/Fields_rzt");
tgrid = h5read("results4.h5","/IRprop/tgrid");
rgrid = h5read("results4.h5","/IRprop/rgrid");
zgrid = h5read("results4.h5","/IRprop/zgrid");

cd(currfold)




kz = 1;
kr = 1024;

Erz1 = squeeze( AllFields(kz,kr,:) );
Erz2 = squeeze( AllFields(kz+1,kr,:) );
figure

plot(tgrid,Erz1);  hold on
plot(tgrid,Erz2)

% plot(tgrid,Erz1,'-*');  hold on
% plot(tgrid,Erz2,'-*')

% plot(Erz1,'-*');  hold on
% plot(Erz2,'-*')

rgrid(kr)
zgrid(kz)


kz = 1;
kt = 1024;

Ezt1 = squeeze( AllFields(kz,:,kt) );
Ezt2 = squeeze( AllFields(kz+1,:,kt) );
Ezt3= squeeze( AllFields(kz+2,:,kt) );
figure

plot(rgrid,Ezt1); hold on
plot(rgrid,Ezt2); hold on
plot(rgrid,Ezt3); hold on

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