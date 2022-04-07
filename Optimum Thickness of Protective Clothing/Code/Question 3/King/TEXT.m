
figure
x=linspace(-0.8,0.8);
y=5*x/2;
[X,Y]=meshgrid(x,y);
Z=20+X.^2-20*cos(2*pi*X)+20+Y.^2-20*cos(2*pi*Y);
mesh(X,Y,Z)
hidden off
hold on
x=G(:,1);y=G(:,2);
z=20+x.^2-20*cos(2*pi*x)+20+y.^2-20*cos(2*pi*y);
scatter3(G(:,1),G(:,2),z,'filled');
hold on
x=G(:,3);y=G(:,4);
z=20+x.^2-20*cos(2*pi*x)+20+y.^2-20*cos(2*pi*y);
scatter3(G(:,3),G(:,4),z,'filled');
hold on
x=G(:,5);y=G(:,6);
z=20+x.^2-20*cos(2*pi*x)+20+y.^2-20*cos(2*pi*y);
scatter3(G(:,5),G(:,6),z,'filled');
hold on
x=G(:,7);y=G(:,8);
z=20+x.^2-20*cos(2*pi*x)+20+y.^2-20*cos(2*pi*y);
scatter3(G(:,7),G(:,8),z,'filled');
hold on
x=G(:,9);y=G(:,10);
z=20+x.^2-20*cos(2*pi*x)+20+y.^2-20*cos(2*pi*y);
scatter3(G(:,9),G(:,10),z,'filled');
hold off