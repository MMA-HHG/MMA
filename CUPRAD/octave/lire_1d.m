function[x,y] = lire_1d(file);

%file
fid= fopen(file,'r');
if fid ~=-1
    buffer=fscanf(fid,'%f');
    fclose(fid);
    [s,dummy]=size(buffer);
    x=zeros(s/2,1);
    y=zeros(s/2,1);
    for n=1:s/2
        x(n,1)=buffer(2*n-1,1);
        y(n,1)=buffer(2*n,1);
    end;
else
    x=zeros(2,1);
    y=zeros(2,1);
end;