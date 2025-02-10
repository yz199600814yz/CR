%%
%确定起终点的像素索引值
%考虑x,y的坐标增量变化,获得初始起点的像素索引值

function s = make_u_v(pointx,pointy)
pointu = fix(pointx);
pointv = fix(pointy);
if pointx - pointu > 0.5 && pointy - pointv <= 0.5
    pointu1 = pointu + 1;
    pointv1 = pointv;
elseif  pointx - pointu > 0.5 && pointy - pointv > 0.5
    pointu1 = pointu + 1;
    pointv1 = pointv + 1;
elseif  pointx - pointu <= 0.5 && pointy - pointv <= 0.5
    pointu1 = pointu;
    pointv1 = pointv;
else
%     pointx - pointu <= 0.5 && pointy - pointv > 0.5;
    pointu1 = pointu;
    pointv1 = pointv + 1;
end
U = pointu1;
V = pointv1;
s =[U,V];
