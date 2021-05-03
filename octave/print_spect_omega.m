function print_spect_omega(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,omegaminhz,omegamaxhz,intmin);

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);   
set(gcf,'color','white','InvertHardCopy','off');

omegaminhz=omegaminhz*1e-15;
omegamaxhz=omegamaxhz*1e-15;

filter = sprintf('%s%s',inpath,'*spect_1d.dat*');
files = dir(filter);
[m,n] = size(files);
k=1;
for n=1:m
    file = sprintf('%s%s',inpath,files(n).name);
    [x,y,u,v,w] = lire_1d_5(file);
    [sx,dummy]=size(x);
    x=x/(tpfs);
    y=y/max(max(y));
    u=u/max(max(u));
    v=v/max(max(v));
    w=unwrap(w);
    graphtitle = strrep(files(n).name,'_spect_1d.dat','');
    z=4*3.1415*n0/(lambdanm*1e-9)*(w0cm/100)^2*str2num(strrep(graphtitle,'_','.'));
    graphtitle = sprintf('%s%0.4g%s','z = ',z,' m');
    subplot(5,2,2*(k-1)+1) %subplot('Position',[0.08 0.85-(k-1)*0.2 0.15 0.12]);
    semilogy(x,y,'Color',[0 0 0]);
    axis ([omegaminhz omegamaxhz intmin 1]);
    if k==5 
       xlabel('\omega [PHz]');
    end;
    ylabel('dE_{Box}/d(\omega)');
    title(graphtitle);
    set(gca,'FontSize',10);
    subplot(5,2,2*(k-1)+2) %subplot('Position',[0.33 0.85-(k-1)*0.2 0.15 0.12]);
    semilogy(x,v,'Color',[0 0 0]);
    axis ([omegaminhz omegamaxhz intmin 1]);
    if k==5 
       xlabel('\omega [PHz]');
    end;
    ylabel('I(r=0,\omega)');title(graphtitle);
    set(gca,'FontSize',10);
    if k==5
        tiffile = sprintf('%s%s%d',outpath,'spect_omega',n/5);
        saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
        print(format,resolution,tiffile);
        close all;
        figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
               'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);   
        set(gcf,'color','white','InvertHardCopy','off');
        k=0;
    end;
    k=k+1;
end;
if k~=1
    tiffile = sprintf('%s%s%d',outpath,'spect_omega',(n+6-k)/5);
    saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
    print(format,resolution,tiffile);
    close all;
end;
