function print_beam(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGw,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);

close all;



 file = sprintf('%s%s',inpath,'fluence.dat');
 [arr,x,y] = lire_2d_bin(file);
 [sx,dummy]=size(x);
 [dummy,sy]=size(y);
 r=zeros(1,sy);
 rmiro=zeros(1,sy);
 for n=1:sy
    for m=1:sx
       if (arr(n,m) > arr(n,1)/2)
           r(1,n)=x(m,1);
       end;
    end;
end;
for n=1:sy
    for m=1:sx
       if (arr(n,m) > arr(n,1)/exp(2))
           rmiro(1,n)=x(m,1);
       end;
    end;
end;
y=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*y;
r=w0cm*10*r; % r_FWHMN
rmiro=w0cm*10*rmiro; %r/e^2

%%% arr(z,r)
% r_stat(ISO)
rgrid=w0cm*10*x;
r_stat=zeros(1,sy);
% r_E % deposed energy
alpha = 0.76;
r_E = zeros(1,sy);
N_r = length(rgrid);

for k1 = 1:sy %over z
	integrand = (arr(k1,:))';
	IntI = cumtrapz(rgrid,integrand);
	IntIM2 = trapz(rgrid,(rgrid.^2).*integrand);
	r_stat(k1) = sqrt(IntIM2/IntI(end));


	% deposed energy
	IntI_search = alpha*IntI(end);

	for k2 = 1:(N_r-1)
		if ( (IntI(k2) >= IntI_search) && (IntI_search < IntI(k2+1)) )
			r_E(k1)=rgrid(k2)+( (IntI_search-IntI(k2))/(IntI(k2+1)-IntI(k2)) )*( rgrid(k2+1)-rgrid(k2) );
			break;
		elseif ( k2 == (N_r-1)) 
			r_E(k1)=-1;
		end
	end
	
end




plot(y,r,'-k');hold on; 
plot(y,rmiro,'-g');hold on;
plot(y,r_stat,'-r'); hold on;
plot(y,r_E,'-b'); hold on;

plot(y,-r,'-k');hold on; 
plot(y,-rmiro,'-g');hold on;
plot(y,-r_stat,'-r'); hold on;
plot(y,-r_E,'-b'); hold on;



%set(gca,'ylim', [-rmaxmm rmaxmm]);

set(gca,'xlim',[zminm zmaxm]);
xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);

% subplot(5,2,4);%subplot('Position',[0.6 0.65 0.35 0.13]);
x=w0cm*10*x;
arr=PcrGw*1.e9/(4*3.1415*(w0cm)^2)*tpfs*1.e-15*arr;
% plot(y,arr(1:sy,1),'Color',[0 0 0]);
% set(gca,'xlim',[zminm zmaxm]);
% xlabel('z [m]');ylabel('FLUENCE(r=0) [J/cm^2]');set(gca,'FontSize',10);

title('Various beam characteristics');
legend({'FWHM','r_{e^2}','r_{rms}','r_{E_{0.76}}'},'Location','best')

saveas(gcf,sprintf('%s%s',outpath,'beam_chars.pdf'));
tiffile = sprintf('%s%s',outpath,'beam_chars');
print(format,resolution,tiffile);
close all;


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

imagesc(yn,xn,arrn);

colorbar;

axis ([zminm zmaxm -rmaxmm rmaxmm],'normal');
xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);
title('FLUENCE [J/cm^2]');


saveas(gcf,sprintf('%s%s',outpath,'beam_fluence.pdf'));
tiffile = sprintf('%s%s',outpath,'beam_fluence');
print(format,resolution,tiffile);
close all;
