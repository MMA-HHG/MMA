function[x,y,z] = lire_1d_3(file);

%file
fid= fopen(file,'r');
if fid ~=-1
buffer=fscanf(fid,'%f');
fclose(fid);
[s,dummy]=size(buffer);
x=zeros(s/3,1);
y=zeros(s/3,1);
z=zeros(s/3,1);
for n=1:s/3
   x(n,1)=buffer(3*n-2,1);
   y(n,1)=buffer(3*n-1,1);
   z(n,1)=buffer(3*n,1);
end;
else
    x=zeros(2,1);
    y=zeros(2,1);
    z=zeros(2,1);
end;