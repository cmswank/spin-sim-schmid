sx=zeros(Bin,Event);
sy=zeros(Bin,Event);
sz=zeros(Bin,Event);
x=zeros(Bin,Event);
y=zeros(Bin,Event);
z=zeros(Bin,Event);
if length(t)>Bin
t=t(2:Bin+1);
end
if length(Sx)>=(Bin)*Event
Sx=Sx(2:end);Sy=Sy(2:end);Sz=Sz(2:end);
X=X(2:end);Y=Y(2:end);Z=Z(2:end);
end
sx(:,1)=Sx(1:Bin);
sy(:,1)=Sy(1:Bin);
sz(:,1)=Sz(1:Bin);
x(:,1)=X(1:Bin);
y(:,1)=Y(1:Bin);
z(:,1)=Z(1:Bin);


for i = 2:Event;
sx(:,i)=Sx(Bin*(i-1)+1:Bin*i);
sy(:,i)=Sy(Bin*(i-1)+1:Bin*i);
sz(:,i)=Sz(Bin*(i-1)+1:Bin*i);
x(:,i)=X(Bin*(i-1)+1:Bin*i);
y(:,i)=Y(Bin*(i-1)+1:Bin*i);
z(:,i)=Z(Bin*(i-1)+1:Bin*i);

end