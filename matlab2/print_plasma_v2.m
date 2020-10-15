function [N_e,N_n] = print_plasma(inpath,outpath,format,resolution,w0cm,tpfs,lambdanm,n0,tminfs,tmaxfs,rmaxmm,rhoccm3,rhontcm3,PrintPart); % N_e is the number of electrons

if PrintPart == 1
	close all;
	figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
	       'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
	set(gcf,'color','white','InvertHardCopy','off');
end;

filter = sprintf('%s%s',inpath,'*_plasma.dat*');
filter
files = dir(filter);
[m,n] = size(files);
k=1;
l=1;
k1=1;



zgrid = zeros(1,m-1);



for n=1:(m-1) % 1:(m-1) 1:10
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

    if n==1 ionmap = zeros((m-1),length(y)); end;
    if n==1 ionmap2 = zeros((m-1),syn); end;
    ionmap(n,:) = (rhoccm3/(rhontcm3*4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2))*arr(:,end); ionmap2(n,:) = arrn(:,end);
    


    % extract z
    graphtitle = strrep(files(n).name,'_plasma.dat','');
    z=4*3.1415*n0/(lambdanm*1e-9)*(w0cm/100)^2*str2num(strrep(graphtitle,'_','.'));
    zgrid(n) = z;
    
    if PrintPart == 1
	    colormap jet;
	   
	    imagesc(xn,yn,100*arrn);...
		
	    colorbar;
	    
		axis ([tminfs tmaxfs -rmaxmm rmaxmm],'normal');...
		if k>=2
		   xlabel('t [fs]');
		end;
		ylabel('r [mm]');...
	    
	   
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
end;
if PrintPart == 1
	if k~=1
	    tiffile = sprintf('%s%s%d',outpath,'plasma',(l+2-k)/2);
	    saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
	    print(format,resolution,tiffile);
	    close all;
	end;
end;


i=find(isnan(ionmap));  ionmap(i)=1;

close all;
figure('Position',[0.393701 0.393701 7.48031 9.84252],'PaperType', 'A4', 'PaperSize',[8.26772,11.69291],...
'PaperPositionMode','manual','PaperPosition',[0 0 8.26772 11.69291]);
set(gcf,'color','white','InvertHardCopy','off');

rgrid = yn;

subplot(2,1,1)

	colormap jet;
   
	imagesc(zgrid,rgrid,100*ionmap2');...
        
	colorbar;

%        axis ([zgrid(1) zgrid(end) -rmaxmm rmaxmm],'normal');...
        xlabel('z [cm]');
        ylabel('r [mm]');...
%        graphtitle = strrep(files(n).name,'_plasma.dat','');
	title('P_{ion}(r,z) [%]');


subplot(2,1,2)
	ionmap2 = rhontcm3 * ionmap2;

	rhocm1=zeros(1,(m-1));
	for k1 = 1:(m-1)		
		grid = rgrid(ceil(syn/2):end);
		polardensity = ionmap2(k1,ceil(syn/2):end)';	
		polardensity = 0.1*grid.*polardensity; % grid in cm
		rhocm1(k1) = 2*3.1415*trapz(0.1*grid,polardensity); % grid in cm, result in cm-1
	end

	plot(zgrid,rhocm1);
	set(gca,'xlim',[zgrid(1) zgrid(end)]);

	xlabel('z [m]');ylabel('\rho [cm^{-1}]');set(gca,'FontSize',10);
	title('Partial electron density');


	N_e = trapz(zgrid,rhocm1);
	N_n = rhontcm3 * 0.1*(zgrid(end)-zgrid(1))*3.1415*(0.1*rgrid(end))^2; % reference number of neutrals in given integrated volume, not very meaningful since r_max is an arbitrary parameter, should be rather related to some fixed cylinder...


%y=rhoccm3/(4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2)*y;
%arrn = (rhoccm3/(rhontcm3*4*3.1415^2*w0cm^2/(lambdanm*1.e-7)^2))*arrn;

%	colormap jet;
%   
%	imagesc(zgrid,rgrid,100*ionmap2);...
%        
%	colorbar;

%        axis ([zgrid(1) zgrid(end) -rmaxmm rmaxmm],'normal');...
%        xlabel('z [cm]');
%        ylabel('r [mm]');...
%        graphtitle = strrep(files(n).name,'_plasma.dat','');




tiffile = sprintf('%s%s',outpath,'ionisationmapv2');
saveas(gcf,sprintf('%s%s',tiffile,'.pdf'));
print(format,resolution,tiffile);
close all;

