function [arr,x,y] = lire_2d_rad_complex(file);

fid= fopen(file,'r','ieee-le');
if fid ~=-1
    fread(fid,1,'int32');entete=fread(fid,2,'int32');fread(fid,1,'int32');
    file
    Nx=entete(1)
    Ny=entete(2)
    arr=complex(zeros(2*Ny-1,Nx),zeros(2*Ny-1,Nx));
    buffer=zeros(1,2*Nx);
    x=zeros(1,Nx);y=zeros(1,2*Ny-1);
    fread(fid,1,'int32');x=fread(fid,Nx,'float32');fread(fid,1,'int32');
    for n=Ny:2*Ny-1
        fread(fid,1,'int32');y(1,n)=fread(fid,1,'float32');fread(fid,1,'int32');
        fread(fid,1,'int32');buffer(1,1:2*Nx)=fread(fid,2*Nx,'float32').';fread(fid,1,'int32');
        arr(n,:)=complex(buffer(1,1:2:2*Nx-1),buffer(1,2:2:2*Nx));
    end;
    for n=1:Ny-1
        y(1,n)=-y(1,2*Ny-n);
        arr(n,:)=arr(2*Ny-n,:);
    end;
    fclose(fid);
    'successfully read'
    clear buffer;

else
    arr=complex(zeros(2,2),zeros(2,2));
    x=zeros(1,2);y=zeros(1,2);
    return;
end;
