function PrintFields(inpath,outpath,w0cm,tpfs,lambdanm,n0,PcrGW,omegauppe,n_resol);

close all;


filter = sprintf('%s%s',inpath,'*field_r_t*');
filter_plasma = sprintf('%s%s',inpath,'*plasma*');
files = dir(filter);
files_plasma = dir(filter_plasma);
[m,n] = size(files);
k=1;
l=1;

zchar_length=6; rchar_length=6;

%cd inpath
%if (exist('fields','dir')==7)
%	rmdir 'fields' s
%end
%mkdir 'fields'


filename = strcat(outpath,'zgrid.dat');
fileID2=fopen(filename,'w')


for n=1:m % 1:m  1:2

	if (n == 1)		
		filename = strcat(outpath,'rgrid.dat');
		fileID1=fopen(filename,'w');
	endif

    file = sprintf('%s%s',inpath,files(n).name);
    [arr,x,y] = lire_2d_rad_complex(file);
    [sy,sx]=size(arr);
    for p=1:sx
        help=exp(-i*omegauppe*x(p));
        arr(:,p)=arr(:,p).*help;
    end;
    
    x=tpfs*x;y=w0cm*10*y;
    arr2=real(sqrt(PcrGW*1e9*3e8*4*3.1415e-7/(4*3.1415*w0cm^2*1e-4*2*n0))*2*arr*1e-9);
    clear arr;
    [sy,sx]=size(arr2);
%    arrplot=arr2(ceil(sy/2),1:sx);



        set(gca,'FontSize',10);
	graphtitle = strrep(files(n).name,'_field_r_t.dat','');
	z=4*3.1415*n0/(lambdanm*1e-9)*(w0cm/100)^2*str2num(strrep(graphtitle,'_','.'));
	fprintf(fileID2,'%i\t%E \n',n,z);

	zchar = sprintf('%i',n);
	k3 = length(zchar);
	if ( k3 < zchar_length)
		for k4 = 1:(zchar_length-k3)
			dumchar(k4)='0';
		end
		zchar=strcat(dumchar,zchar);
		clear dumchar;
	endif

	k2 = 0;	k5 = 0;
	for k1 = ceil(sy/2):sy %ceil(sy/2):sy  ceil(sy/2):(ceil(sy/2)+2)
		k2 = k2+1;
		if ( mod(k2-1,n_resol) == 0 )
			k5 = k5+1;
			rchar = sprintf('%i',k5);
			k3 = length(rchar);
			if ( k3 < rchar_length)
				for k4 = 1:(rchar_length-k3)
					dumchar(k4)='0';
				end
				rchar=strcat(dumchar,rchar);
				clear dumchar;
			endif

			filename = strcat(outpath,'field_z_',zchar,'_r_',rchar,'.dat');
			fileID = fopen(filename,'w');
			fx = arr2(k1,1:sx);
			for k4 = 1:length(x)
			fprintf(fileID,'%E \n',fx(k4));
			end
			fclose(fileID);
			if (n == 1)		
				fprintf(fileID1,'%i\t%E \n',k5,y(k1));
				if (k1 == sy )
					fclose(fileID1);
				endif
				if (k1 == ceil(sy/2) )
					filename = strcat(outpath,'tgrid.dat');
					fileID = fopen(filename,'w');
					for k4 = 1:length(x)
						fprintf(fileID,'%E \n',x(k4));
					end
					fclose(fileID);
				endif
			endif
		endif
	end
end;

fclose(fileID2);
