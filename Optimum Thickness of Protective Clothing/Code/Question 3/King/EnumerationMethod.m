clear,clc
l4=5.5;
ET=65;time=60;
l2=(0.6:0.1:25);
m=length(l2);
EvaluationMatrix=zeros(m,3);
fprintf('评价矩阵开始记录\n');
for i=1:m
    fprintf('记录到第%d代\n',i);
    [FinalTemperature,DurationTime]=Q2_CrankNicolson(l2(i),l4,ET,time);
    EvaluationMatrix(i,1)=l2(i);
    EvaluationMatrix(i,2)=FinalTemperature;
    EvaluationMatrix(i,3)=DurationTime;
end
fprintf('评价矩阵记录结束\n'); 

MeanEva=mean(EvaluationMatrix);
StandardizationMatrix=zeros(m,3);
StandardizationMatrix(:,1)=(EvaluationMatrix(:,1)-MeanEva(1))./(sqrt(sum((EvaluationMatrix(:,1)-MeanEva(1)).^2)/(m-1)));
StandardizationMatrix(:,2)=(EvaluationMatrix(:,2)-MeanEva(2))./(sqrt(sum((EvaluationMatrix(:,2)-MeanEva(2)).^2)/(m-1)));
StandardizationMatrix(:,3)=(EvaluationMatrix(:,3)-MeanEva(3))./(sqrt(sum((EvaluationMatrix(:,3)-MeanEva(3)).^2)/(m-1)));

MeanSta=mean(StandardizationMatrix);
delta(1)=sqrt(sum((StandardizationMatrix(:,1)-MeanSta(1)).^2)/(m-1));
delta(2)=sqrt(sum((StandardizationMatrix(:,2)-MeanSta(2)).^2)/(m-1));
delta(3)=sqrt(sum((StandardizationMatrix(:,3)-MeanSta(3)).^2)/(m-1));

V=zeros(1,3);
V(1)=delta(1)/MeanSta(1);
V(2)=delta(2)/MeanSta(2);
V(3)=delta(3)/MeanSta(3);

W=zeros(3,3);
for i=1:3
    W(i,i)=V(i)/sum(V);
end

PJ=StandardizationMatrix*W;

SolutionBest=[min(PJ(:,1)),max(PJ(:,2)),max(PJ(:,3))];
SolutionWorst=[max(PJ(:,1)),min(PJ(:,2)),min(PJ(:,3))];

dbest=zeros(m,1);
dworst=zeros(m,1);
for i=1:m
    dbest(i)=sqrt((PJ(i,1)-SolutionBest(1))^2+(PJ(i,2)-SolutionBest(2))^2+(PJ(i,3)-SolutionBest(3))^2);
    dworst(i)=sqrt((PJ(i,1)-SolutionWorst(1))^2+(PJ(i,2)-SolutionWorst(2))^2+(PJ(i,3)-SolutionWorst(3))^2);
end

C=dworst./(dworst+dbest);

for i=1:m
    if EvaluationMatrix(i,1)>25||EvaluationMatrix(i,1)<0.6
        C(i)=0;
    end
    if EvaluationMatrix(i,2)>=47
        C(i)=0;
    end
    if EvaluationMatrix(i,3)>=300
        C(i)=0;
    end
end

[~,I]=max(C);
BestAnswer=EvaluationMatrix(I,:);
BestLength=BestAnswer(1)
BestFinalTemperature=BestAnswer(2)
BestDurationTime=BestAnswer(3)

figure
plot(1:3600,UcnM(220,:),'LineWidth',2)
grid on
xlabel('时间 / 秒','fontsize',15)
ylabel('温度 / 摄氏度','fontsize',15)
hold on
plot([3600-52,3600-52],[36,45],'--r','LineWidth',1.5)
gtext('开始超过44摄氏度--->')


figure
plot(1:256,UcnM(:,3600),'LineWidth',2)
grid on
xlabel('深度 / 毫米','fontsize',15)
ylabel('温度 / 摄氏度','fontsize',15)
hold on
scatter(6,UcnM(6,3600),'filled');
hold on
scatter(129,UcnM(129,3600),'filled');
hold on
scatter(165,UcnM(165,3600),'filled');
hold on
scatter(220,UcnM(220,3600),'filled');
hold on
scatter(256,UcnM(256,3600),'filled');
hold on
x=[6,6];y=[35,75]; 
plot(x,y,'--k','LineWidth',1.5);
hold on
x=[129,129];y=[35,75]; 
plot(x,y,'--k','LineWidth',1.5);
hold on
x=[165,165];y=[35,75]; 
plot(x,y,'--k','LineWidth',1.5);
x=[220,220];y=[35,75]; 
plot(x,y,'--k','LineWidth',1.5);
text(2,50,'I层');
text(80,50,'II层');
text(140,50,'III层');
text(190,50,'IV层');
text(235,50,'皮肤内侧');

