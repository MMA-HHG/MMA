function [omegauppe,w0cm,tpfs,lambdanm,n0,PcrGW,rhoccm3,rhontcm3,rhontabscm3,sigmak,KK] = lire_listing(file);

sigmak=0;
KK=1;
file
fid= fopen(file,'r');
if fid ~=-1
  line=fgetl(fid);
  while line ~=-1 
      line=fgetl(fid);
  end;
  line=fgetl(fid);
  while line ~=-1
      if line(1:9)=='beamwaist'
         w0cm=sscanf(line,'%*s%f')
      end;
      if line(1:14)=='pulse duration'
         tpfs=sscanf(line,'%*s%*s%f')
      end;     
      if line(1:10)=='wavelenght'
         lambdanm=1.e7*sscanf(line,'%*s%*s%f')
      end;
      if line(1:16)=='refractive index'
         n0=sscanf(line,'%*s%*s%*s%*s%f')
      end;
      if line(1:16)=='number of photon'
         KK=sscanf(line,'%*s%*s%*s%*s%*s%*s%*s%f')
      end;
      if line(1:19)=='coefficient for MPI'
         sigmak=sscanf(line,'%*s%*s%*s%f')
      end;
      if line(1:15)=='critical plasma'
         rhoccm3=sscanf(line,'%*s%*s%*s%f')
      end;
      if line(1:17)=='effective density'
         if line(38)=='s'
            rhontcm3=sscanf(line,'%*s%*s%*s%*s%*s%f')
	 end;
      end;
      if line(1:17)=='Density of absorb'
         rhontabscm3=sscanf(line,'%*s%*s%*s%*s%f')'
      end;
      line=fgetl(fid);
  end;
  line=fgetl(fid);
  while line ~=-1
      if line(1:14)=='critical power'
         PcrGW=1.e-9*sscanf(line,'%*s%*s%f')'
      end;
      line=fgetl(fid);
  end;
  line=fgetl(fid);
  line=fgetl(fid);
  while line ~=-1
      if line(1:42)=='adimensionned central frequency of the box'
         omegauppe=sscanf(line,'%*s%*s%*s%*s%*s%*s%f')'
      end;
      line=fgetl(fid);
  end;
  
  fclose(fid);
else 
  close all;
return;
end;
