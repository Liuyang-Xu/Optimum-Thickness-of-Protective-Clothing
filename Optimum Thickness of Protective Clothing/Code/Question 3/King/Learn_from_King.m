function [VG] = Learn_from_King(Gbest,King,radius,King_number,human_number,dimension,city_number,List,c3,c4,c5)

%LEARN_FROM_KING �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
distance=repmat(King,1,city_number)-Gbest;                                          %��ǰ��������֮��ľ���
a=1:city_number;
a(a==King_number)=[];
Vstack=zeros(1,dimension*city_number);
VG=zeros(human_number,dimension*city_number);
%����������ѧϰ������
for i=a
    %��һά��
    if abs(distance(1,dimension*i-1))>radius(1,dimension*i-1)
        Vstack(1,dimension*i-1)=c3(i)*distance(1,dimension*i-1).*rand;
    end
    
    %if abs(distance(1,dimension*i-1))<c3(1,dimension*i-1)
        %Vstack(1,dimension*i-1)=-0.5*(distance(1,dimension*i-1)).*rand;
    %end
    
    %�ڶ�ά��
    if abs(distance(1,dimension*i))>radius(1,dimension*i)
        Vstack(1,dimension*i)=c3(i)*distance(1,dimension*i).*rand;
    end
    
    %if abs(distance(1,dimension*i))<c3(1,dimension*i)
        %Vstack(1,dimension*i)=-0.5*(distance(1,dimension*i)).*rand;
    %end
    
    if (abs(distance(1,dimension*i))<radius(1,dimension*i)&&abs(distance(1,dimension*i-1))<radius(1,dimension*i-1))
        Vstack(1,dimension*i)=c4(i)*(distance(1,dimension*i)-radius(1,dimension*i)).*rand;
        Vstack(1,dimension*i-1)=c4(i)*(distance(1,dimension*i-1)-radius(1,dimension*i-1)).*rand;
    end
    
    
    
    VG(List(1,i),dimension*i-1:dimension*i)=Vstack(1,dimension*i-1:dimension*i);
end
VG(List(1,King_number),dimension*King_number-1:dimension*King_number)=c5.*rand(1,2);


end

