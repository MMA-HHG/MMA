function print_profile(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGw,rhoccm3,rhontcm3,sigmak,KK,zminm,zmaxm,rmaxmm,tminfs,tmaxfs,intmintwcm2,intmaxtwcm2,omegaminhz,omegamaxhz);

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);   
set(gcf,'color','white','InvertHardCopy','off');



file = sprintf('%s%s',inpath,'intensity_onax_t.dat');
[arr,x,y] = lire_2d_bin(file);
[sx,dummy]=size(x);
[dummy,sy]=size(y);
y=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*y;
x=tpfs*x;
arr=PcrGw/(1000*4*3.1415*w0cm^2)*arr.*arr;

subplot(2,1,1);%subplot('Position',[0.1 0.84 0.35 0.13]);
file = sprintf('%s%s',inpath,'peakmax.dat');
[xalt,yalt] = lire_1d(file);
xalt= (4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*xalt;
yalt=PcrGw/(1000*4*3.1415*w0cm^2)*yalt.*yalt;
[xp,yp]=reduce(xalt,yalt);
plot(xp,yp,'Color',[0 0 0]);hold on;
plot(y,arr(1:sy,sx/2),'Color',[0 0 0],'LineStyle','--');
axis([zminm zmaxm -inf inf]);
xlabel('z [m]');ylabel('I [TW/cm^2]');set(gca,'FontSize',10);




subplot(2,1,2);

file = sprintf('%s%s',inpath,'rhomax.dat');
[x,y] = lire_1d(file);
x=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*x;
% y=rhoccm3/(4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)*y;
y=rhoccm3/(rhontcm3*4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)*y;
[xp,yp]=reduce(x,y);

% semilogy(xp,yp,'Color',[0 0 0]);

plot(xp,100*yp);

set(gca,'xlim',[zminm zmaxm]);
xlabel('z [m]');ylabel('Ionisation [%]');set(gca,'FontSize',10);


saveas(gcf,sprintf('%s%s',outpath,'intensity.pdf'));
tiffile = sprintf('%s%s',outpath,'intensity');
print(format,resolution,tiffile);





% figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
%        'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);   
% set(gcf,'color','white','InvertHardCopy','off');
% 
%  subplot(5,2,2);%subplot('Position',[0.6 0.84 0.35 0.13]);
%  file = sprintf('%s%s',inpath,'fluence.dat');
%  [arr,x,y] = lire_2d_bin(file);
%  [sx,dummy]=size(x);
%  [dummy,sy]=size(y);
%  r=zeros(1,sy);
%  rmiro=zeros(1,sy);
%  for n=1:sy
%     for m=1:sx
%        if (arr(n,m) > arr(n,1)/2)
%            r(1,n)=x(m,1);
%        end;
%     end;
% end;
% for n=1:sy
%     for m=1:sx
%        if (arr(n,m) > arr(n,1)/exp(2))
%            rmiro(1,n)=x(m,1);
%        end;
%     end;
% end;
% y=(4.*3.1415*n0*(w0cm/100)^2/(lambdanm*1.e-9))*y;
% r=w0cm*10*r;
% rmiro=w0cm*10*rmiro;
% plot(y,r,'Color',[0 0 0]);hold on;plot(y,-r,'Color',[0 0 0]);hold on;
% plot(y,rmiro,'Color',[0 0 0],'LineStyle','--');hold on;
% plot(y,-rmiro,'Color',[0 0 0],'LineStyle','--');
% set(gca,'xlim',[zminm zmaxm]);
% xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);
% 
% 
% 
% 
% 
% 
% subplot(5,2,4);%subplot('Position',[0.6 0.46 0.35 0.13]);
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
% 
% colormap jet;
% 
% imagesc(yn,xn,arrn);
% axis ([zminm zmaxm -rmaxmm rmaxmm],'normal');
% xlabel('z [m]');ylabel('r [mm]');set(gca,'FontSize',10);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% saveas(gcf,sprintf('%s%s',outpath,'profile.pdf'));
% tiffile = sprintf('%s%s',outpath,'profile');
% print(format,resolution,tiffile);
% % close all;