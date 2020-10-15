%function post_process;

close all; clear all;

path=pwd
inpath=sprintf('%s%s',path,'/')
outpath=sprintf('%s%s',path,'/figures2/')

% cd ../octave

format='-djpeg'
resolution='-r300'
zminm=0.
zmaxm=.05
rmaxmm=.3
tminfs=-100
tmaxfs=100
% omegaminhz=1.8e15
% omegamaxhz=3e15
lambdaminnm=500
lambdamaxnm=1100
intmin=1.e-6
intminlin=1.e-6
intmintwcm2=0
intmaxtwcm2=100

omegaminhz=2*3.1415*3e17/lambdamaxnm
omegamaxhz=2*3.1415*3e17/lambdaminnm

PrintPart = 0
printvecz = [1:10:71, 81:2:101]
Nplot_endplane = 80


[omegauppe,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,rhontabscm3,sigmak,KK] = lire_listing(sprintf('%s%s',inpath,'listing'));


print_profile(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
print_beam_separ(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
print_beam(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
print_test(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz); % review
print_spect_lambda(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intminlin,lambdaminnm,lambdamaxnm); % spectrum

print_slice_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,tminfs,tmaxfs,rmaxmm,printvecz); % intensity

print_slice_field_endplane(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,tminfs,tmaxfs,rmaxmm,Nplot_endplane);  %intensity endplane

print_slice_field_real(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,omegauppe,tminfs,tmaxfs,rmaxmm); % field
%[N_e2,N_n2]=print_plasma_v2(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm,rhoccm3,rhontcm3,PrintPart);
[N_e,N_n]=print_plasma(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm,rhoccm3,rhontcm3,zmaxm,zminm);


print_spect_lambdalog(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin,lambdaminnm,lambdamaxnm);

print_spect_omega(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin,printvecz);



%print_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegauppe,tminfs,tmaxfs,rmaxmm);

display('Number of electrons')
%[N_e N_n]
%[N_e2 N_n2]



%print_spect_omega_endplane(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin,Nplot_endplane);
