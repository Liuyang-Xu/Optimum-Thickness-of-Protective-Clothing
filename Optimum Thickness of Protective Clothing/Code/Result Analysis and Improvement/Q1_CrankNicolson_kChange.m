close all
clear
clc

LeftBoundry=xlsread('HumanSkinTemperature.xlsx','b3:b5403');

l1=0.6;
l2=6;
l3=3.6;
l4=5;

p=[300,862,74.2,1.18]*10^(-9);%�ܶ�
C=[1377,2100,1726,1005];%����
K=[0.082,0.37,0.045,0.028]*10^(-3);%�ȴ�����

l=l1+l2+l3+l4;%��λmm

h=0.1;
k=0.01;

time=90;%��λmin
time=time*60;%��λs

Mx=ceil(l/h)+1;%������x���ϵĽڵ����������0��ĩβ��
Nt=ceil(time/k)+1;%������t���ϵĽڵ����������0��ĩβ��

alpha=zeros(1,4);
for i=1:4
    alpha(i)=K(i)/(C(i)*p(i));
end

r=zeros(1,4);
for i=1:4
    r(i)=alpha(i)*k/(2*(h*h));%�����
end

%================================================================

%׷�Ϸ����
%An*U(n+1)=Bn*Un+en
 
%��An
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

AnU=zeros(1,Mx-3);%�ϴζԽ�
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

AnL=zeros(1,Mx-3);%�´ζԽ�
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

%��Bn
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

BnU=zeros(1,Mx-3);%�ϴζԽ�
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

BnL=zeros(1,Mx-3);%�´ζԽ�
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

%׷�Ϸ�
U0=zeros(1,Mx);%��ֵ
for x=1:Mx
    U0(x)=Q1_InitialConditions_kChange(x);
end
 
UcnM=zeros(Mx-2,Nt-1);
 
%��U1  An*U(n+1)=Bn*Un+en
en=zeros(Mx-2,1);
en(1)=r(1)*Q1_BoundaryConditions_kChange(0,1,LeftBoundry)+r(1)*Q1_BoundaryConditions_kChange(0,2,LeftBoundry);
en(Mx-2)=r(4)*Q1_BoundaryConditions_kChange(1,1,LeftBoundry)+r(4)*Q1_BoundaryConditions_kChange(1,2,LeftBoundry);

UcnM(:,1)=Cn*U0(2:Mx-1)'+An\en;

%ǰ50s
for t=2:5000
    n=t-1;
    fprintf('%d\n',n);
    en=zeros(Mx-2,1);
    en(1)=r(1)*Q1_BoundaryConditions_kChange(0,n,LeftBoundry)+r(1)*Q1_BoundaryConditions_kChange(0,(n+1),LeftBoundry);
    en(Mx-2)=r(4)*Q1_BoundaryConditions_kChange(1,n,LeftBoundry)+r(4)*Q1_BoundaryConditions_kChange(1,(n+1),LeftBoundry);
    UcnM(:,n+1)=Cn*UcnM(:,n)+An\en;
end

UcnMfirst=UcnM(:,1:5000);

%======================================================



%======================================================

k=1;

Resttime=time-50;%��λs

newNt=ceil(Resttime/k)+1;%������t���ϵĽڵ����������0��ĩβ��

newU0=UcnM(:,5000);%��ֵ

newUcnM=zeros(Mx-2,newNt-1);

%��U1  An*U(n+1)=Bn*Un+en
newen=zeros(Mx-2,1);
newen(1)=r(1)*75+r(1)*75;
newen(Mx-2)=r(4)*LeftBoundry(50)+r(4)*LeftBoundry(51);

newUcnM(:,1)=Cn*newU0+An\newen;

%����ʱ����
for t=2:newNt-1
    n=t-1;
    fprintf('newis%d\n',n);
    newen=zeros(Mx-2,1);
    newen(1)=r(1)*75+r(1)*75;
    newen(Mx-2)=r(4)*LeftBoundry(50+n)+r(4)*LeftBoundry(50+(n+1));
    newUcnM(:,n+1)=Cn*newUcnM(:,n)+An\newen;
end


%================================================================

xgrid = h:h:l-h;
tgrid = k/100:k/100:50;
[X, T] = meshgrid(tgrid,xgrid);
figure
surf(X,T,UcnMfirst)
xlabel('ʱ�� / ��','fontsize',15)
ylabel('��� / ����','fontsize',15)
zlabel('�¶� / ���϶�','fontsize',15)
shading interp

hold on

xgrid = h:h:l-h;
tgrid = 50+k:k:time;
[X, T] = meshgrid(tgrid,xgrid);
surf(X,T,newUcnM)
shading interp

%======================================================
