%function post_process;

inpath=pwd
inpath=sprintf('%s%s',inpath,'/')
outpath=sprintf('%s%s',inpath,'figures/')

cd ~/CUPRAD/octave

format='-djpeg'
resolution='-r300'
zminm=0.
zmaxm=2.
rmaxmm=.3
tminfs=-100
tmaxfs=100
omegaminhz=0e15
omegamaxhz=1e16
intmin=1.e-6
intmintwcm2=0
intmaxtwcm2=100

[omegauppe,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,rhontabscm3,sigmak,KK] = lire_listing(sprintf('%s%s',inpath,'listing'));

print_iso(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);
print_spect_omega(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin);
print_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegauppe,tminfs,tmaxfs,rmaxmm);
print_plasma(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm);
print_slice_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,tminfs,tmaxfs,rmaxmm);
print_slice_field_real(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,omegauppe,tminfs,tmaxfs,rmaxmm);
