function [N_e,N_n] = print_plasma(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm,rhoccm3,rhontcm3,zmaxm,zminm); % N_e is the number of electrons


close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
set(gcf,'color','white','InvertHardCopy','off');

subplot(2,1,1)

file = sprintf('%s%s',inpath,'plasmachannel.dat');
[arr,x,y] = lire_2d_bin(file);
[sx,dummy]=size(x);
[dummy,sy]=size(y);
y=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*y;
x=w0cm*10*x;
arr=rhoccm3/(rhontcm3*4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)*arr;
sxn=999;
syn=1000;
xn=zeros(sxn,1);
for n=1:floor(sxn/2)
    xn(floor(sxn/2)+1+n,1)=n*rmaxmm/floor(sxn/2);
end;
for n=1:floor(sxn/2)
    xn(n,1)=-xn(sxn+1-n,1);
end;
yn=zeros(1,syn);
for n=1:syn
    yn(1,n)=(n-1)*(zmaxm-zminm)/(syn-1)+zminm;
end;
arrn=zeros(sxn,syn);
arrn(floor(sxn/2)+1:sxn,1:syn)=interp2(y,x,arr',yn,xn(floor(sxn/2)+1:sxn,1));
for n=1:floor(sxn/2)
    arrn(n,1:syn)=arrn(sxn+1-n,1:syn);
end;
i=find(isnan(arrn));
arrn(i)=0; %Usual contrast = 1
colormap jet;
imagesc(yn,xn,100*arrn);
axis ([zminm zmaxm -rmaxmm rmaxmm],'normal');
xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);




	
   
%	imagesc(zgrid,rgrid,100*ionmap2');...
        
	colorbar;

%        axis ([zgrid(1) zgrid(end) -rmaxmm rmaxmm],'normal');...
        xlabel('z [cm]');
        ylabel('r [mm]');...
%        graphtitle = strrep(files(n).name,'_plasma.dat','');
	title('P_{ion}(r,z) [%]');


subplot(2,1,2)
	m = length(yn);
	ionmap2 = rhontcm3 * arrn';
	rgrid = xn; zgrid = yn;

	rhocm1=zeros(1,m);
	for k1 = 1:m	
		grid = rgrid(ceil(sxn/2):end);
		polardensity = ionmap2(k1,ceil(sxn/2):end)';	
		polardensity = 0.1*grid.*polardensity; % grid in cm
		rhocm1(k1) = 2*3.1415*trapz(0.1*grid,polardensity); % grid in cm, result in cm-1
	end

	plot(zgrid,rhocm1);
	set(gca,'xlim',[zgrid(1) zgrid(end)]);

	xlabel('z [m]');ylabel('\rho [cm^{-1}]');set(gca,'FontSize',10);
	title('Partial electron density');


	N_e = trapz(zgrid,rhocm1);
	N_n = rhontcm3 * 0.1*(zgrid(end)-zgrid(1))*3.1415*(0.1*rgrid(end))^2; % reference number of neutrals in given integrated volume, not very meaningful since r_max is an arbitrary parameter, should be rather related to some fixed cylinder...




tiffile = sprintf('%s%s',outpath,'ionisationmap');
saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
print(format,resolution,tiffile);
close all;

