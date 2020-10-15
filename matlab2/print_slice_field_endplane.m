function print_slice_field(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,tminfs,tmaxfs,rmaxmm,Nplot);

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
set(gcf,'color','white','InvertHardCopy','off');

filter = sprintf('%s%s',inpath,'*field_r_t*');
filter_plasma = sprintf('%s%s',inpath,'*plasma*');
files = dir(filter);
files_plasma = dir(filter_plasma);
[m,n] = size(files);
k=1;
l=1;

n = m;
file = sprintf('%s%s',inpath,files(n).name);
[arr,x,y] = lire_2d_rad_complex(file);
arr2=abs(arr);
clear arr;
arr=PcrGW/(1000*4*3.1415*w0cm^2)*arr2.*arr2;
clear arr2;

[sy,sx]=size(arr);
dumint = floor(ceil(sy/2)/(Nplot+1))
vecplot = ones(1,Nplot);
for k1 = 1: Nplot
	vecplot(k1) = ceil(sy/2) + (k1-1) * dumint;
end

%vecplot = [2048 2049 2070 2080 2090 2100 2420 2606 2792 3350]

x=tpfs*x;y=w0cm*10*y;




for k1=1:Nplot

%display('loop')

    subplot(5,2,(k-1)+1) %subplot('Position',[0.08 0.85-(k-1)*0.2 0.15 0.12]);
    
    

    
    arrplot=arr(vecplot(k1),1:sx);
    plot(x,arrplot,'Color',[0 0 0]);
    set(gca,'xlim',[tminfs tmaxfs]);
    if k==5
       xlabel('t [fs]');
    end;
    

    graphtitle = strrep(files(n).name,'_field_r_t.dat','');
    z=4*3.1415*n0/(lambdanm*1e-9)*(w0cm/100)^2*str2num(strrep(graphtitle,'_','.'));

    graphtitle = sprintf('%s%0.2g%s','z = ',100*z,' cm');...
    ylabel(strcat('I(',graphtitle,') [TW/cm^2]'));
    graphtitle = sprintf('%s%3.0f%s','r = ',1e3*y(vecplot(k1)),' {\mu}m');...

        title(graphtitle);...
        set(gca,'FontSize',10);

%	if ismember(n,plotvecz)
%			filename = strcat(outpath,strcat('intensity_r0_',num2str(n),'.dat');
%			fileID = fopen(filename,'w');
%			for k4 = 1: (length(x)-1)
%			fprintf(fileID,'%E \t %E \n',x(k4), arrplot(k4));
%			end
%			fprintf(fileID,'%E \t %E ',x(length(x)), arrplot(length(x)) );
%			fclose(fileID);
%	elseif n == m
%			filename = strcat(outpath,'intensity_r0_end.dat');
%			fileID = fopen(filename,'w');
%			for k4 = 1: (length(x)-1)
%			fprintf(fileID,'%E \t %E \n',x(k4), arrplot(k4));
%			end
%			fprintf(fileID,'%E \t %E ',x(length(x)), arrplot(length(x)) );
%			fclose(fileID);
%	end


%    subplot(5,2,2*(k-1)+2) %subplot('Position',[0.33 0.85-(k-1)*0.2 0.15 0.12]);
%    arrplot=arr(1:sy,ceil(sx/2));
%    plot(y,arrplot,'Color',[0 0 0]);
%    set(gca,'xlim',[-rmaxmm rmaxmm]);
%    if k==5
%       xlabel('r [mm]');
%    end;
%    ylabel('I(t=0) [TW/cm^2]');
%    title(graphtitle);
%    clear arr;
%    clear arrplot;


    if k==10
        tiffile = sprintf('%s%s%d',outpath,'field_slices_endplane',l/10);
        saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
        print(format,resolution,tiffile);
        close all;
        figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
        set(gcf,'color','white','InvertHardCopy','off');
        k=0;
    end;
    k=k+1;
    l=l+1;
end;
if k~=1
    tiffile = sprintf('%s%s%d',outpath,'field_slices_endplane',(l+5-k)/5);
    saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
    print(format,resolution,tiffile);
    close all;
end;


filename = strcat(outpath,strcat('intensities_endplane.dat'));
fileID = fopen(filename,'w');
for k1=1:sx
	for k2 = 1:(Nplot-1)
		fprintf(fileID,'%E \t', arr(vecplot(k2),k1) );	
	end
	fprintf(fileID,'%E \n', arr(vecplot(Nplot),k1) );
end
fclose(fileID);

filename = strcat(outpath,strcat('tgrid_endplane.dat'));
fileID = fopen(filename,'w');
for k1=1:sx
	fprintf(fileID,'%E \n', x(k1) );	
end
fclose(fileID);

filename = strcat(outpath,strcat('rgrid_endplane.dat'));
fileID = fopen(filename,'w');
for k1 = 1:(Nplot-1)
	fprintf(fileID,'%E \n', y(vecplot(k1)) );	
end
fclose(fileID);



