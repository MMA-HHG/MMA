function [xp,yp]=reduce(x,y)

xp=x;
yp=y;
[sx dummy]=size(xp);
n=2;
while n<=sx
    if (xp(n-1,1)>=xp(n,1))
        xp(n,:)=[];
        yp(n,:)=[];
        n=n-1;
    end;
    [sx dummy]=size(xp);
    n=n+1;
end;