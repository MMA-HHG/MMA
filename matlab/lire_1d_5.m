function[x,y,z,v,w] = lire_1d_5(file);

%file
fid= fopen(file,'r');
if fid ~=-1
buffer=fscanf(fid,'%f');
fclose(fid);
[s,dummy]=size(buffer);
x=zeros(s/5,1);
y=zeros(s/5,1);
z=zeros(s/5,1);
v=zeros(s/5,1);
w=zeros(s/5,1);
for n=1:s/5
   x(n,1)=buffer(5*n-4,1);
   y(n,1)=buffer(5*n-3,1);
   z(n,1)=buffer(5*n-2,1);
   v(n,1)=buffer(5*n-1,1);
   w(n,1)=buffer(5*n,1);
end;
else
    x=zeros(2,1);
    y=zeros(2,1);
    z=zeros(2,1);
    v=zeros(2,1);
    w=zeros(2,1);
end;