function [arr,x,y] = lire_2d_bin(file);

fid= fopen(file,'r','ieee-le');
if fid ~=-1
  fread(fid,1,'int32');entete=fread(fid,2,'int32');fread(fid,1,'int32');
  %Nx=entete(1);Ny=entete(2);  
  Nx=entete(1);Ny=entete(2);
  arr=zeros(Ny,Nx);
  x=zeros(1,Nx);y=zeros(1,Ny);
  fread(fid,1,'int32');x=fread(fid,Nx,'float32');fread(fid,1,'int32');
  for n=1:Ny
    fread(fid,1,'int32');y(1,n)=fread(fid,1,'float32');fread(fid,1,'int32');
    fread(fid,1,'int32');arr(n,:)=fread(fid,Nx,'float32')';fread(fid,1,'int32');
  end;
  fclose(fid);
else
  arr=zeros(2,2);
  x=zeros(2,1);y=zeros(1,2);
  x(2,1)=1;
  y(1,2)=1;
return;
end;
