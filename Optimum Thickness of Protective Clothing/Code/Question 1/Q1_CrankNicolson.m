close all
clear
clc

LeftBoundry=xlsread('HumanSkinTemperature.xlsx','b3:b5403');

l1=0.6;
l2=6;
l3=3.6;
l4=5;

p=[300,862,74.2,1.18]*10^(-9);%密度
C=[1377,2100,1726,1005];%比热
K=[0.082,0.37,0.045,0.028]*10^(-3);%热传导率

l=l1+l2+l3+l4;%单位mm

h=0.1;
k=1;

time=90;%单位min
time=time*60;%单位s

Mx=ceil(l/h)+1;%网格在x轴上的节点个数（算上0和末尾）
Nt=ceil(time/k)+1;%网格在t轴上的节点个数（算上0和末尾）

alpha=zeros(1,4);
for i=1:4
    alpha(i)=K(i)/(C(i)*p(i));
end

r=zeros(1,4);
for i=1:4
    r(i)=alpha(i)*k/(2*(h*h));%网格比
end

%================================================================

%追赶法求解
%An*U(n+1)=Bn*Un+en
 
%求An
An=zeros(1,Mx-2);
for i=1:ceil((l1/h)-1)
    An(i)=1+2*r(1);
end
for i=ceil((l1/h)+1):ceil((l1+l2)/h-1)
    An(i)=1+2*r(2);
end
for i=ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1)
    An(i)=1+2*r(3);
end
for i=ceil((l1+l2+l3)/h+1):ceil(l/h-1)
    An(i)=1+2*r(4);
end

AnU=zeros(1,Mx-3);%上次对角
for i=1:ceil((l1/h)-1)
    AnU(i)=-r(1);
end
for i=ceil((l1/h)+1):ceil((l1+l2)/h-1)
    AnU(i)=-r(2);
end
for i=ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1)
    AnU(i)=-r(3);
end
for i=ceil((l1+l2+l3)/h+1):ceil(l/h-1)
    AnU(i)=-r(4);
end

AnL=zeros(1,Mx-3);%下次对角
for i=2:ceil((l1/h)-1)
    AnL(i)=-r(1);
end
for i=ceil((l1/h)+1):ceil((l1+l2)/h-1)
    AnL(i)=-r(2);
end
for i=ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1)
    AnL(i)=-r(3);
end
for i=ceil((l1+l2+l3)/h+1):ceil(l/h-1)
    AnL(i)=-r(4);
end
An=diag(An,0)+diag(AnU(1:end-1),1)+diag(AnL(2:end),-1);

for i=ceil(l1/h)
    An(i,i-2)=K(1)/2;
    An(i,i-1)=-2*K(1);
    An(i,i)=3*(K(1)+K(2))/2;
    An(i,i+1)=-2*K(2);
    An(i,i+2)=K(2)/2;
end
for i=ceil((l1+l2)/h)
    An(i,i-2)=K(2)/2;
    An(i,i-1)=-2*K(2);
    An(i,i)=3*(K(2)+K(3))/2;
    An(i,i+1)=-2*K(3);
    An(i,i+2)=K(3)/2;
end
for i=ceil((l1+l2+l3)/h)
    An(i,i-2)=K(3)/2;
    An(i,i-1)=-2*K(3);
    An(i,i)=3*(K(3)+K(4))/2;
    An(i,i+1)=-2*K(4);
    An(i,i+2)=K(4)/2;
end

%================================================================

%求Bn
Bn=zeros(1,Mx-2);
for i=1:ceil((l1/h)-1)
    Bn(i)=1-2*r(1);
end
for i=ceil((l1/h)+1):ceil((l1+l2)/h-1)
    Bn(i)=1-2*r(2);
end
for i=ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1)
    Bn(i)=1-2*r(3);
end
for i=ceil((l1+l2+l3)/h+1):ceil(l/h-1)
    Bn(i)=1-2*r(4);
end

BnU=zeros(1,Mx-3);%上次对角
for i=1:ceil((l1/h)-1)
    BnU(i)=r(1);
end
for i=ceil((l1/h)+1):ceil((l1+l2)/h-1)
    BnU(i)=r(2);
end
for i=ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1)
    BnU(i)=r(3);
end
for i=ceil((l1+l2+l3)/h+1):ceil(l/h-1)
    BnU(i)=r(4);
end

BnL=zeros(1,Mx-3);%下次对角
for i=2:ceil((l1/h)-1)
    BnL(i)=r(1);
end
for i=ceil((l1/h)+1):ceil((l1+l2)/h-1)
    BnL(i)=r(2);
end
for i=ceil((l1+l2)/h+1):ceil((l1+l2+l3)/h-1)
    BnL(i)=r(3);
end
for i=ceil((l1+l2+l3)/h+1):ceil(l/h-1)
    BnL(i)=r(4);
end
Bn=diag(Bn,0)+diag(BnU(1:end-1),1)+diag(BnL(2:end),-1);

for i=ceil(l1/h)
    for j=(i-2):(i+2)
        Bn(i,j)=0;
    end
end
for i=ceil((l1+l2)/h)
    for j=(i-2):(i+2)
        Bn(i,j)=0;
    end
end
for i=ceil((l1+l2+l3)/h)
    for j=(i-2):(i+2)
        Bn(i,j)=0;
    end
end

%================================================================

Cn=An\Bn;

%================================================================

%追赶法
U0=zeros(1,Mx);%初值
for x=1:Mx
    U0(x)=Q1_InitialConditions(x);
end
 
UcnM=zeros(Mx-2,Nt-1);
 
%求U1  An*U(n+1)=Bn*Un+en
en=zeros(Mx-2,1);
en(1)=r(1)*Q1_BoundaryConditions(0,1,LeftBoundry)+r(1)*Q1_BoundaryConditions(0,2,LeftBoundry);
en(Mx-2)=r(4)*Q1_BoundaryConditions(1,1,LeftBoundry)+r(4)*Q1_BoundaryConditions(1,2,LeftBoundry);

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
    fprintf('%d\n',n);
    en=zeros(Mx-2,1);
    en(1)=r(1)*Q1_BoundaryConditions(0,n,LeftBoundry)+r(1)*Q1_BoundaryConditions(0,(n+1),LeftBoundry);
    en(Mx-2)=r(4)*Q1_BoundaryConditions(1,n,LeftBoundry)+r(4)*Q1_BoundaryConditions(1,(n+1),LeftBoundry);
    UcnM(:,n+1)=Cn*UcnM(:,n)+An\en;
end


U=zeros(Mx,Nt);
U(:,1)=U0';
U(1,:)=ones(1,Nt)*75;
U(Mx,:)=LeftBoundry;
for i=2:Mx-1
    for j=2:Nt
        U(i,j)=UcnM(i-1,j-1);
    end
end

%================================================================

xgrid = 0:h:l;
tgrid = 0:k:time;
[X, T] = meshgrid(tgrid,xgrid);
figure
surf(X,T,U)
xlabel('Time / s','fontsize',15)
ylabel('Depth/ mm','fontsize',15)
zlabel('Temperature / ℃','fontsize',15)
shading interp

%ax^2+bx
%[a1,b1,a2,b2,x2]=solve('a1*(102^2)+b1*102=65.1203','a1*(152^2)+b1*152=48.0800','2*152*a2+b2=(2*a1*152+b1)*2.8000e-05/(0.035*10^(-3))','a2*(152^2)+b2*152=48.0800','a2*(x2^2)+b2*x2=37')

%a/x+b
%[a1,b1,a2,b2,x2]=solve('a1/102+b1=65.1203','a1/152+b1=48.0800','2.8000e-05-a1/(102^2)+a2/(152^2)=0.035*10^(-3)','a2/152+b2=48.0800','a2/x2+b2=37')

%alogx+b=====================
%[a1,b1,a2,b2,x2]=solve('a1*log(102)+b1=65.1023','a1*log(152)+b1=48.0800','2.8000e-05*a1/102=0.035*10^(-3)*a2/152','a2*log(152)+b2=48.0800','a2*log(x2)+b2=37')

%ax^2+c
%[a1,b1,a2,b2,x2]=solve('a1*(102^2)+b1=65.1203','a1*(152^2)+b1=48.0800','2*152*a2=(2*a1*152)*2.8000e-05/(0.035*10^(-3))','a2*(152^2)+b2=48.0800','a2*(x2^2)+b2=37')

%ax^3+bx^2
%[a1,b1,a2,b2,x2]=solve('a1*(102^3)+b1*(102^2)=65.1203','a1*(152^3)+b1*(152^2)=48.0800','2.8000e-05*(3*a1*(102^2)+2*b1*102)=0.035*10^(-3)*(3*a2*(152^2)+2*b2*152)','a2*(152^3)+b2*(152^2)=48.0800','a2*(x2^3)+b2*(x2^2)=37')

%ax^3+bx
%[a1,b1,a2,b2,x2]=solve('a1*(102^3)+b1*(102)=65.1203','a1*(152^3)+b1*(152)=48.0800','2.8000e-05*(3*a1*(102^2)+b1)=0.035*10^(-3)*(3*a2*(152^2)+b2)','a2*(152^3)+b2*(152)=48.0800','a2*(x2^3)+b2*(x2)=37')

%ae^bx
%[a1,b1,a2,b2,x2]=solve('a1*exp(b1*102)=65.1203','a1*exp(b1*102)=48.0800','2.8000e-05*a1*b1*exp(102)=0.2*10^(-3)*a2*b2*exp(152)','a2*exp(b2*152)=48.0800','a2*exp(b2*x2)=37')

figure
plot(1:153,U(:,5401),'LineWidth',2)
grid on
xlabel('Depth/ mm','fontsize',15)
ylabel('Temperature / ℃','fontsize',15)
hold on
scatter(7,U(7,5401),'filled');
hold on
scatter(67,U(67,5401),'filled');
hold on
scatter(103,U(103,5401),'filled');
hold on
scatter(153,U(153,5401),'filled');
hold on
x=[7,7];y=[45,75]; 
plot(x,y,'--k','LineWidth',1.5);
hold on
x=[67,67];y=[45,75]; 
plot(x,y,'--k','LineWidth',1.5);
hold on
x=[103,103];y=[45,75]; 
plot(x,y,'--k','LineWidth',1.5);
text(2,50,'Layer I');
text(30,50,'Layer II');
text(80,50,'Layer III');
text(130,50,'Layer IV');

figure
plot(1:5401,U(3,:),'LineWidth',2);
hold on
plot(1:5401,U(35,:),'LineWidth',2);
hold on
plot(1:5401,U(80,:),'LineWidth',2);
hold on
plot(1:5401,U(130,:),'LineWidth',2);
grid on
xlabel('Time / s','fontsize',15)
ylabel('Temperature / ℃','fontsize',15)
legend('Temperature Curve of Layer I ','Temperature Curve of Layer II ','Temperature Curve of Layer III ','Temperature Curve of Layer IV ')







