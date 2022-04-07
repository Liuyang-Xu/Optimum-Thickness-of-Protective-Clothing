function [FinalTemperature,DurationTime,UcnM]=Q2_CrankNicolson(l2,l4,ET,time)

l1=0.6;
%l2;
l3=3.6;
%l4;
%l5=76.2;
%l5=105.4;
l5=3.7;

p=[300,862,74.2,1.18,2300]*10^(-9);%密度
C=[1377,2100,1726,1005,350];%比热
K=[0.082,0.37,0.045,0.028,0.035]*10^(-3);%热传导率

l=l1+l2+l3+l4+l5;%单位mm

h=0.1;
k=1;%EnumerationMethod

%time;%单位min
time=time*60;%转化为单位s

Mx=ceil(l/h)+1;%网格在x轴上的节点个数（算上0和末尾）
Nt=ceil(time/k)+1;%网格在t轴上的节点个数（算上0和末尾）

alpha=K./(C.*p);
r=alpha*k/(2*(h*h));%网格比

%================================================================

%追赶法求解
%An*U(n+1)=Bn*Un+en

%求An
An=zeros(1,Mx-2);
An(1:ceil((l1/h)-1))=1+2*r(1);
An(ceil((l1/h)+1):ceil((l1+l2)/h-1))=1+2*r(2);
An(ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1))=1+2*r(3);
An(ceil((l1+l2+l3)/h+1):ceil((l1+l2+l3+l4)/h-1))=1+2*r(4);
An(ceil((l1+l2+l3+l4)/h+1):ceil(l/h-1))=1+2*r(5);

%上次对角
AnU=zeros(1,Mx-3);
AnU(1:ceil((l1/h)-1))=-r(1);
AnU(ceil((l1/h)+1):ceil((l1+l2)/h-1))=-r(2);
AnU(ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1))=-r(3);
AnU(ceil((l1+l2+l3)/h+1):ceil((l1+l2+l3+l4)/h-1))=-r(4);
AnU(ceil((l1+l2+l3+l4)/h+1):ceil(l/h-1))=-r(5);

%下次对角
AnL=zeros(1,Mx-3);
AnL(2:ceil((l1/h)-1))=-r(1);
AnL(ceil((l1/h)+1):ceil((l1+l2)/h-1))=-r(2);
AnL(ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1))=-r(3);
AnL(ceil((l1+l2+l3)/h+1):ceil((l1+l2+l3+l4)/h-1))=-r(4);
AnL(ceil((l1+l2+l3+l4)/h+1):ceil(l/h-1))=-r(5);
An=diag(An,0)+diag(AnU(1:end-1),1)+diag(AnL(2:end),-1);

i=ceil(l1/h);
An(i,i-2)=K(1)/2;
An(i,i-1)=-2*K(1);
An(i,i)=3*(K(1)+K(2))/2;
An(i,i+1)=-2*K(2);
An(i,i+2)=K(2)/2;
i=ceil((l1+l2)/h);
An(i,i-2)=K(2)/2;
An(i,i-1)=-2*K(2);
An(i,i)=3*(K(2)+K(3))/2;
An(i,i+1)=-2*K(3);
An(i,i+2)=K(3)/2;
i=ceil((l1+l2+l3)/h);
An(i,i-2)=K(3)/2;
An(i,i-1)=-2*K(3);
An(i,i)=3*(K(3)+K(4))/2;
An(i,i+1)=-2*K(4);
An(i,i+2)=K(4)/2;
i=ceil((l1+l2+l3+l4)/h);
An(i,i-2)=K(4)/2;
An(i,i-1)=-2*K(4);
An(i,i)=3*(K(4)+K(5))/2;
An(i,i+1)=-2*K(5);
An(i,i+2)=K(5)/2;

%================================================================

%求Bn
Bn=zeros(1,Mx-2);
Bn(1:ceil((l1/h)-1))=1-2*r(1);
Bn(ceil((l1/h)+1):ceil((l1+l2)/h-1))=1-2*r(2);
Bn(ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1))=1-2*r(3);
Bn(ceil((l1+l2+l3)/h+1):ceil((l1+l2+l3+l4)/h-1))=1-2*r(4);
Bn(ceil((l1+l2+l3+l4)/h+1):ceil(l/h-1))=1-2*r(5);

%上次对角
BnU=zeros(1,Mx-3);
BnU(1:ceil((l1/h)-1))=r(1);
BnU(ceil((l1/h)+1):ceil((l1+l2)/h-1))=r(2);
BnU(ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1))=r(3);
BnU(ceil((l1+l2+l3)/h+1):ceil((l1+l2+l3+l4)/h-1))=r(4);
BnU(ceil((l1+l2+l3+l4)/h+1):ceil(l/h-1))=r(5);

%下次对角
BnL=zeros(1,Mx-3);
BnL(2:ceil((l1/h)-1))=r(1);
BnL(ceil((l1/h)+1):ceil((l1+l2)/h-1))=r(2);
BnL(ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1))=r(3);
BnL(ceil((l1+l2+l3)/h+1):ceil((l1+l2+l3+l4)/h-1))=r(4);
BnL(ceil((l1+l2+l3+l4)/h+1):ceil(l/h-1))=r(5);
Bn=diag(Bn,0)+diag(BnU(1:end-1),1)+diag(BnL(2:end),-1);

i=ceil(l1/h);
Bn(i,(i-2):(i+2))=0;
i=ceil((l1+l2)/h);
Bn(i,(i-2):(i+2))=0;
i=ceil((l1+l2+l3)/h);
Bn(i,(i-2):(i+2))=0;
i=ceil((l1+l2+l3+l4)/h);
Bn(i,(i-2):(i+2))=0;

%================================================================

Cn=An\Bn;

%================================================================

%追赶法
U0=zeros(1,Mx);%初值
U0(1)=ET;
U0(2:Mx)=37;

UcnM=zeros(Mx-2,Nt-1);

%求U1  An*U(n+1)=Bn*Un+en
en=zeros(Mx-2,1);
en(1)=r(1)*ET+r(1)*ET;
en(Mx-2)=r(5)*37+r(5)*37;

UcnM(:,1)=Cn*U0(2:Mx-1)'+An\en;

%{
bn=Bn*U0(2:Mx-1)';
for i = 2:Mx-2
    An(i,i-1) = An(i,i-1)/An(i-1,i-1);
    An(i,i) = An(i,i) - An(i-1,i) * An(i,i-1);
    bn(i) = bn(i) - bn(i-1) * An(i,i-1);
end
 
     x(Mx-2) = bn(Mx-2) / An(i,i);
 
for i =Mx-2-1:-1:1
    x(i) = (bn(i) - An(i,i+1) * x(i+1)) / An(i,i);
end
x=x';
UcnM(:,1)=x+en;
%}

for t=2:Nt-1
    n=t-1;
    en=zeros(Mx-2,1);
    en(1)=r(1)*ET+r(1)*ET;
    en(Mx-2)=r(5)*37+r(5)*37;
    UcnM(:,n+1)=Cn*UcnM(:,n)+An\en;
end

%================================================================

FinalTemperature=UcnM(ceil((l1+l2+l3+l4)/h),Nt-1);
%MaxTemperature=max(UcnM(ceil((l1+l2+l3+l4)/h),:));
for FirstTime=1:Nt-1
    if UcnM(ceil((l1+l2+l3+l4)/h),FirstTime)>44
        break
    end
end
DurationTime=((Nt-1)-FirstTime)*k;

%================================================================

%{
xgrid = h:h:l-h;
tgrid = k:k:time;
[X,T] = meshgrid(tgrid,xgrid);
figure
surf(X,T,UcnM)
xlabel('Time','fontsize',15)
ylabel('Length','fontsize',15)
zlabel('Temperature','fontsize',15)
shading interp
%}


