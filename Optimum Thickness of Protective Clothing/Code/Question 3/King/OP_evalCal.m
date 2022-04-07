function [EvalValue] = OP_evalCal(G,human_number,dimension,city_number)
%OP_EVALCAL 此处显示有关此函数的摘要
%   此处显示详细说明
%   输出结果是列向量

EvalValue=zeros(human_number,city_number);               %初始化

ET=80;
time=30;
Ideal_Point_Worst=[1,1,1,1];
Ideal_Point_Best=[0.6/25,0.6/6.4,37/80,0];

G=ceil(G*10)*0.1;

for i=1:city_number
    for j=1:human_number
        l2=G(j,dimension*i-1);
        l4=G(j,dimension*i);
        [FinalTemperature,DurationTime]=Q2_CrankNicolson(l2,l4,ET,time);
        if DurationTime>300||FinalTemperature>47||FinalTemperature<37||l2>25||l2<0.6||l4>6.4||l4<0.6
            EvalValue(j,i)=0;
        else
            dbest=sqrt( (Ideal_Point_Best(1)-l2/25)^2+ (Ideal_Point_Best(2)-l4/6.4)^2 + (Ideal_Point_Best(3)-FinalTemperature/80)^2  + (Ideal_Point_Best(4)-DurationTime/1800)^2   );
            dworst=sqrt( (Ideal_Point_Worst(1)-l2/25)^2+ (Ideal_Point_Worst(2)-l4/6.4)^2 + (Ideal_Point_Worst(3)-FinalTemperature/80)^2  + (Ideal_Point_Worst(4)-DurationTime/1800)^2   );

            EvalValue(j,i)=dworst/(dbest+dworst);
        end
    end
end

end

