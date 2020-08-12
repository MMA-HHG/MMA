%function post_process;

close all; clear all;

path=pwd
inpath=sprintf('%s%s',path,'/')
outpath=sprintf('%s%s',path,'/figures2/')

% cd ../octave

format='-djpeg'
resolution='-r300'
zminm=0.
zmaxm=.2
rmaxmm=.3
tminfs=-100
tmaxfs=100
% omegaminhz=1.8e15
% omegamaxhz=3e15
lambdaminnm=700
lambdamaxnm=900
intmin=1.e-6
intmintwcm2=0
intmaxtwcm2=100

omegaminhz=2*3.1415*3e17/lambdamaxnm
omegamaxhz=2*3.1415*3e17/lambdaminnm

[omegauppe,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,rhontabscm3,sigmak,KK] = lire_listing(sprintf('%s%s',inpath,'listing'));

print_profile(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
print_beam(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
% print_test(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
% % print_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegauppe,tminfs,tmaxfs,rmaxmm);
print_plasma(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm,rhoccm3,rhontcm3);
print_spect_lambda(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin,lambdaminnm,lambdamaxnm);
print_spect_omega(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin);
% % print_slice_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,tminfs,tmaxfs,rmaxmm);
print_slice_field_real(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,omegauppe,tminfs,tmaxfs,rmaxmm);
