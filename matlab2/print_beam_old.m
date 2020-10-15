function print_beam(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGw,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);   
set(gcf,'color','white','InvertHardCopy','off');



 subplot(2,1,1);%subplot('Position',[0.6 0.84 0.35 0.13]);
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
for k1 = 1:sy %over z
	integrand = (arr(k1,:))';
	IntI = trapz(rgrid,integrand);
	IntIM2 = trapz(rgrid,(rgrid.^2).*integrand);
	r_stat(k1) = sqrt(IntIM2/IntI);	
end



plot(y,r,'Color',[0 0 0]);hold on;plot(y,-r,'Color',[0 0 0]);hold on;
plot(y,rmiro,'Color',[0 0 0],'LineStyle','--');hold on;
plot(y,-rmiro,'Color',[0 0 0],'LineStyle','--');

plot(y,r_stat,'--r');

set(gca,'xlim',[zminm zmaxm]);
xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);

% subplot(5,2,4);%subplot('Position',[0.6 0.65 0.35 0.13]);
x=w0cm*10*x;
arr=PcrGw*1.e9/(4*3.1415*(w0cm)^2)*tpfs*1.e-15*arr;
% plot(y,arr(1:sy,1),'Color',[0 0 0]);
% set(gca,'xlim',[zminm zmaxm]);
% xlabel('z [m]');ylabel('FLUENCE(r=0) [J/cm^2]');set(gca,'FontSize',10);

subplot(2,1,2);%subplot('Position',[0.6 0.46 0.35 0.13]);
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

% disp('test1')
% subplot(5,2,3);%subplot('Position',[0.1 0.65 0.35 0.13]);
% disp('trig2')
% 
% file = sprintf('%s%s',inpath,'energy.dat')
% [x,y,z] = lire_1d_3(file);
% x=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*x;
% y=PcrGw*1.e9/(4*3.1415)*tpfs*1.e-15*1.e6*y;
% z=PcrGw*1.e9/(4*3.1415)*tpfs*1.e-15*1.e6*z;
% [xp,yp]=reduce(x,y);
% plot(xp,yp,'Color',[0 0 0]);hold on;
% [xp,zp]=reduce(x,z);
% plot(xp,zp,'Color',[0 0 0],'LineStyle','--');
% set(gca,'xlim',[zminm zmaxm]);
% xlabel('z [m]');ylabel('ENERGY [\muJ]');set(gca,'FontSize',10);
% 
% subplot(5,2,5);%subplot('Position',[0.1 0.46 0.35 0.13]);
% file = sprintf('%s%s',inpath,'rhomax.dat');
% [x,y] = lire_1d(file);
% x=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*x;
% y=rhoccm3/(4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)*y;
% [xp,yp]=reduce(x,y);
% semilogy(xp,yp,'Color',[0 0 0]);
% set(gca,'xlim',[zminm zmaxm]);
% xlabel('z [m]');ylabel('\rho [cm^{-3}]');set(gca,'FontSize',10);
% 
% subplot(5,2,7);%subplot('Position',[0.1 0.27 0.35 0.13]);
% file = sprintf('%s%s',inpath,'powmax.dat');
% [x,y] = lire_1d(file);
% x=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*x;
% y=y/(4*3.1415);
% [xp,yp]=reduce(x,y);
% plot(xp,yp,'Color',[0 0 0]);
% set(gca,'xlim',[zminm zmaxm]);
% xlabel('z [m]');ylabel('P/P_{cr}');set(gca,'FontSize',10);
% 
% subplot(5,2,9);%subplot('Position',[0.1 0.08 0.35 0.13]);
% file = sprintf('%s%s',inpath,'ionisation_table.dat');
% [x,y,z] = lire_1d_3(file);
% [s,dummy]=size(x);
% if (s == 2)
%     x = intmintwcm2:(intmaxtwcm2-intmintwcm2)/100:intmaxtwcm2;
%     y = sigmak*(1.e12*x).^KK;
% else
%     x=PcrGw/(1000*4*3.1415*w0cm^2)*x;
%     y=rhoccm3/(4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)/(rhontcm3*tpfs*1e-15)*y;
% end;
% plot(x,y,'Color',[0 0 0]);
% axis([intmintwcm2 intmaxtwcm2 0 1.1*y(min(floor(intmaxtwcm2*s/x(s)),s))]);
% xlabel('I [TW/cm^2]');ylabel('IONISATION RATE [s^{-1}]');set(gca,'FontSize',10);
% 
% %fid= fopen('O2_PPT.dat', 'wt');
% %[s,dummy]=size(x);
% %for n=1:s
% %    fprintf(fid,'%e\t%e\n',x(n),y(n)*1.e-15);
% %end;
% %fclose(fid);
%    
% file = sprintf('%s%s',inpath,'intensity_onax_t.dat');
% [arr,x,y] = lire_2d_bin(file);
% [sx,dummy]=size(x);
% [dummy,sy]=size(y);
% y=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*y;
% x=tpfs*x;
% arr=PcrGw/(1000*4*3.1415*w0cm^2)*arr.*arr;
% 
% subplot(5,2,1);%subplot('Position',[0.1 0.84 0.35 0.13]);
% file = sprintf('%s%s',inpath,'peakmax.dat');
% [xalt,yalt] = lire_1d(file);
% xalt= (4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*xalt;
% yalt=PcrGw/(1000*4*3.1415*w0cm^2)*yalt.*yalt;
% [xp,yp]=reduce(xalt,yalt);
% plot(xp,yp,'Color',[0 0 0]);hold on;
% plot(y,arr(1:sy,sx/2),'Color',[0 0 0],'LineStyle','--');
% axis([zminm zmaxm -inf inf]);
% xlabel('z [m]');ylabel('I [TW/cm^2]');set(gca,'FontSize',10);
% 
% subplot(5,2,8);%subplot('Position',[0.6 0.27 0.35 0.13]);
% sxn=1000;
% syn=2000;
% xn=zeros(sxn,1);
% for n=1:sxn
%     xn(n,1)=(n-1)*(tmaxfs-tminfs)/(sxn-1)+tminfs;
% end;
% yn=zeros(1,syn);
% for n=1:syn
%     yn(1,n)=(n-1)*(zmaxm-zminm)/(syn-1)+zminm;
% end;
% arrn=zeros(sxn,syn);
% %y2=y(1:2:sy);%if arr too large
% %arr2=arr(1:2:sy,:);
% %arrn=interp2(y2,x,arr2',yn,xn);
% arrn=interp2(y,x,arr',yn,xn);
% i=find(isnan(arrn));
% arrn(i)=0; %Usual contrast = 1
% imagesc(yn,xn,arrn);
% axis ([zminm zmaxm tminfs tmaxfs],'normal');
% xlabel('z [m]');ylabel('t [fs]');set(gca,'FontSize',10);
% 
% %omegaminhz=omegaminhz*1e-15;
% %omegamaxhz=omegamaxhz*1e-15;
% %subplot(5,2,10);%subplot('Position',[0.6 0.08 0.3 0.13]);
% %file = sprintf('%s%s',inpath,'nomega.dat');
% %[x,y,z] = lire_1d_3(file);
% %x=x/(tpfs);
% %[s,dummy]=size(x);
% %[haxes,hline1,hline2] = plotyy(x,y,x,z,'plot','semilogy');
% %axes(haxes(1));
% %set(gca,'xlim',[omegaminhz omegamaxhz]);xlabel('\omega [PHz]');set(gca,'FontSize',10,'YColor',[0 0 0]);
% %set(gca,'ylim',[y(floor(1+(max(omegaminhz,x(1))-x(1))*s/(x(s)-x(1)))) y(floor((min(omegamaxhz,x(s))-x(1))*s/(x(s)-x(1))))]);
% %ylabel('n');
% %axes(haxes(2));
% %ylabel('\kappa');
% %set(gca,'xlim',[omegaminhz omegamaxhz]);xlabel('\omega [PHz]');set(gca,'FontSize',10,'YColor',[0 0 0]);
% %set(hline1,'Color',[0 0 0]);
% %set(hline2,'LineStyle','--','Color',[0 0 0]);
% 
% subplot(5,2,10);%subplot('Position',[0.6 0.84 0.35 0.13]);
% file = sprintf('%s%s',inpath,'plasmachannel.dat');
% [arr,x,y] = lire_2d_bin(file);
% [sx,dummy]=size(x);
% [dummy,sy]=size(y);
% y=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*y;
% x=w0cm*10*x;
% arr=rhoccm3/(4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)*arr;
% sxn=999;
% syn=1000;
% xn=zeros(sxn,1);
% for n=1:floor(sxn/2)
%     xn(floor(sxn/2)+1+n,1)=n*rmaxmm/floor(sxn/2);
% end;
% for n=1:floor(sxn/2)
%     xn(n,1)=-xn(sxn+1-n,1);
% end;
% yn=zeros(1,syn);
% for n=1:syn
%     yn(1,n)=(n-1)*(zmaxm-zminm)/(syn-1)+zminm;
% end;
% arrn=zeros(sxn,syn);
% arrn(floor(sxn/2)+1:sxn,1:syn)=interp2(y,x,arr',yn,xn(floor(sxn/2)+1:sxn,1));
% for n=1:floor(sxn/2)
%     arrn(n,1:syn)=arrn(sxn+1-n,1:syn);
% end;
% i=find(isnan(arrn));
% arrn(i)=0; %Usual contrast = 1
% imagesc(yn,xn,arrn);
% axis ([zminm zmaxm -rmaxmm rmaxmm],'normal');
% xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);

saveas(gcf,sprintf('%s%s',outpath,'beam.pdf'));
tiffile = sprintf('%s%s',outpath,'beam');
print(format,resolution,tiffile);
close all;
