function print_plasma(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm,rhoccm3,rhontcm3);

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
set(gcf,'color','white','InvertHardCopy','off');

filter = sprintf('%s%s',inpath,'*_plasma.dat*');
filter
files = dir(filter);
[m,n] = size(files);
k=1;
l=1;
k1=1;
for n=1:m
    file = sprintf('%s%s',inpath,files(n).name);
    [arr,x,y] = lire_2d_rad_real(file);
    subplot(2,1,k); %subplot('Position',[0.08+(k-1-4*floor((k-1)/4))*0.25 0.85-floor((k-1)/4)*0.2 0.15 0.12]);
    x=tpfs*x;y=w0cm*10*y;
    sxn=600;
    syn=299;
    xn=zeros(1,sxn);
    for p=1:sxn
        xn(1,p)=(p-1)*(tmaxfs-tminfs)/(sxn-1)+tminfs;
    end;
    yn=zeros(syn,1);
    for p=1:floor(syn/2)
        yn(floor(syn/2)+1+p,1)=p*rmaxmm/floor(syn/2);
    end;
    for p=1:floor(syn/2)
        yn(p,1)=-yn(syn+1-p,1);
    end;
    arrn=zeros(syn,sxn);
    arrn(floor(syn/2)+1:syn,1:sxn)=interp2(x,y,arr,xn,yn(floor(syn/2)+1:syn,1));
    for p=1:floor(syn/2)
        arrn(p,1:sxn)=arrn(syn+1-p,1:sxn);
    end;
    
    arrn = (rhoccm3/(rhontcm3*4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2))*arrn;
    i=find(isnan(arrn));
    arrn(i)=1;
    
    
    colormap jet;
   
    imagesc(xn,yn,100*arrn);...
        
    colorbar;
    
        axis ([tminfs tmaxfs -rmaxmm rmaxmm],'normal');...
        if k>=2
           xlabel('t [fs]');
        end;
        ylabel('r [mm]');...
        graphtitle = strrep(files(n).name,'_plasma.dat','');
    z=4*3.1415*n0/(lambdanm*1e-9)*(w0cm/100)^2*str2num(strrep(graphtitle,'_','.'));
    graphtitle = sprintf('%s%0.4g%s','z = ',z,' m');...
        title(graphtitle);...
        set(gca,'FontSize',10);
    if k==2
        tiffile = sprintf('%s%s%d',outpath,'plasma',l/2);
        saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
        print(format,resolution,tiffile);
        close all;
        figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
        set(gcf,'color','white','InvertHardCopy','off');
        k=0
    end;
    k=k+1;
    l=l+1;
end;
if k~=1
    tiffile = sprintf('%s%s%d',outpath,'plasma',(l+2-k)/2);
    saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
    print(format,resolution,tiffile);
    close all;
end;
