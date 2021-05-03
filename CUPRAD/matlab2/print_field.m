function print_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegauppe,tminfs,tmaxfs,rmaxmm);

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
set(gcf,'color','white','InvertHardCopy','off');

filter = sprintf('%s%s',inpath,'*field_r_t*');
files = dir(filter);
[m,n] = size(files);
k=1;
l=1;
for n=1:m
    file = sprintf('%s%s',inpath,files(n).name);
    [arr,x,y] = lire_2d_rad_complex(file);%lire_2d_rad_complex(file);
    [sy,sx]=size(arr);
    for p=1:sx
        help=exp(-i*omegauppe*x(p));
        arr(:,p)=arr(:,p).*help;
    end;
    subplot(5,2,k); %subplot('Position',[0.08+(k-1-4*floor((k-1)/4))*0.25 0.85-floor((k-1)/4)*0.2 0.15 0.12]);
    x=tpfs*x;y=w0cm*10*y;
    tlo=1;
    tup=sx;
    for p=1:sx
        if x(p)<=tminfs
            tlo=p;
        end;
        if x(p)>=tmaxfs
            tup=p;
            break;
        end;
    end;
    xlo=1;
    xup=sy;
    for p=1:sy
        if y(p)<=-rmaxmm
            xlo=p;
        end;
        if y(p)>=rmaxmm
            xup=p;
            break;
        end;
    end;
    arr2=abs(arr(xlo:xup,tlo:tup));
    colormap jet;
    imagesc(x(tlo:tup),y(xlo:xup),arr2);...
        axis ([tminfs tmaxfs -rmaxmm rmaxmm]);...
        if k>=9
           xlabel('t [fs]');
        end;
        ylabel('r [mm]');...
        graphtitle = strrep(files(n).name,'_field_r_t.dat','');...
        z=4*3.1415*n0/(lambdanm*1e-9)*(w0cm/100)^2*str2num(strrep(graphtitle,'_','.'));...
        graphtitle = sprintf('%s%0.4g%s','z = ',z,' m');...
        title(graphtitle);...
        set(gca,'FontSize',10);
    if k==10
        tiffile = sprintf('%s%s%d',outpath,'field',l/10);
        saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
        print(format,resolution,tiffile);
        close all;
        figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
               'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
        set(gcf,'color','white','InvertHardCopy','off');
        k=0;
    end;
    clear arr arr2;
    k=k+1;
    l=l+1;
end;
if k~=1
    tiffile = sprintf('%s%s%d',outpath,'field',(l+10-k)/10);
    saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
    print(format,resolution,tiffile);
    close all;
end;
