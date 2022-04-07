function [G,V] = people_learn(G,Pbest,Gbest,V,w,c1,c2,c3,c4,c5,range,Vmax,Vmin,human_number,dimension,city_number,radius,King_number,King,List)
%people_learn �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
%   ��Ⱥѧϰ
Vmax=repmat(Vmax,human_number,1);
Vmin=repmat(Vmin,human_number,1);
Group=repmat(Gbest,human_number,1);
c1=repmat(kron(c1,ones(1,dimension)),human_number,1);
c2=repmat(kron(c2,ones(1,dimension)),human_number,1);

w=repmat(kron(w,ones(1,dimension)),human_number,1);
V=w.*V+c1.*rand(human_number,dimension*city_number).*(Pbest-G)+c2.*rand(human_number,dimension*city_number).*(Group-G);     %�ٶȸ���

VG=Learn_from_King(Gbest,King,radius,King_number,human_number,dimension,city_number,List,c3,c4,c5);
V=V+VG;

for i=1:human_number
    for j=1:dimension*city_number
        if abs(V(i,j))>abs(Vmax(i,j))
            V(i,j)=sign(V(i,j))*Vmax(i,j);
        end
        if abs(V(i,j))<abs(Vmin(i,j))
            V(i,j)=sign(V(i,j))*Vmin(i,j);
        end
    end
end
%�������С�ٶ������ٶ�����
V(List(1,King_number),dimension*King_number-1:dimension*King_number)=[0,0];%����ֹ����ǰ
G=G+V;                                                      %λ�ø���
K=zeros(human_number,dimension*city_number,2);
K(:,:,1)=repmat(range(:,2)',human_number,city_number);
K(:,:,2)=repmat(range(:,1)',human_number,city_number);
for i=1:human_number
    for j=1:dimension*city_number
        if G(i,j)>K(i,j,1)
            G(i,j)=K(i,j,1);
        else
            if G(i,j)<K(i,j,2)
                G(i,j)=K(i,j,2);
            end
        end
    end
end
end
%Խ�����
