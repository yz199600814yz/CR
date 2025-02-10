clc;
clear;
close all;
%% 数据准备，影像和控制点读取
% im1 = uint8(imread('影像\ZW-3M-5624.tiff'));
im1 = imread('影像\SS-1m-2170.tiff');
GCP = load('控制点\SS-1m-2170.txt');
% GCP = GCP(:,2:3);
pointnum = size(GCP,1);
%% 获取控制点矩形窗口
% step1 将点坐标转换为像素坐标并取邻域内的灰度值和坐标索引
win = 1;               % 窗口大小阈值
ww = cell(pointnum,1); % 创建空数组用于存储窗口索引的灰度值
Index1 = cell(pointnum,1); %创建空数组用于存储窗口索引值
uv = [];                   %控制点整数坐标
for i = 1:pointnum 
    %temp = make_u_v(GCP(i,1),GCP(i,2));
    temp = make_u_v(GCP(i,2),GCP(i,3));
    tempw = im1(temp(:,2)-win:temp(:,2)+win,temp(:,1)-win:temp(:,1)+win);
    % 取所取控制点的矩形区域索引
    tempx = [temp(:,1)-win:temp(:,1)+win];
    tempy = [temp(:,2)-win:temp(:,2)+win];
    index = [];
    for j = 1:size(tempx,2)
        aa = [];
        for k = 1:size(tempy,2)
            tempm = [tempx(:,j) tempy(:,k)];
            aa = [aa;tempm];
        end
        index = [index,aa];
    end
    num = size(tempx,2);   
    bb = zeros(1,num);
    for ii = 1:size(bb,2)
        bb(ii) = 1;
    end
    dim1 = bb;
    dim2 = bb.*2;
    Index = mat2cell(index,dim1,dim2);
    ww{i,:} = tempw;   % 矩形框灰度元胞数组
    Index1{i,:} = Index;  %矩形框像素索引元胞数组
    uv = [uv;temp];
end
%%  对窗口内像素进行亚像素分割
% step2 取像素中心点并设定步长对像素进行分割，通过权重函数对分割像素快进行赋值
np = 10;              %设定分割数
length = 1/np;                 %分割步长
% 对矩阵窗口内的像素进行亚像素分割
% 1.取元胞数组的左上角和右下角索
SUBindex = cell(pointnum,1);
for i = 1:size(Index1,1)
    winpt = Index1{i,:};             %3*3像素索引元胞数组
    [m,n] = size(winpt);
    minuv = winpt{1,1};
    maxuv = winpt{m,n};
    x0 = [minuv(:,1):length:maxuv(:,1)];
    y0 = [minuv(:,2):length:maxuv(:,2)];
    tempX = [];
    for j = 1:size(x0,2)
        tempY = [];
        for k = 1:size(y0,2)
            % 进行亚像素细分
            temp_sub = [x0(:,j),y0(:,k)];
            tempY = [tempY;temp_sub];
        end
        tempX = [tempX tempY];
    end
    subnum = size(tempX ,1);   
    cc = zeros(1,subnum);
    for ii = 1:size(cc,2)
        cc(ii) = 1;
    end
    dim3 = cc;
    dim4 = cc.*2;
    subIndex = mat2cell(tempX,dim3,dim4);
    SUBindex{i,:} = subIndex;
end
%% 2. 通过权重函数计算每个点的灰度值
[suby,FKG,FKID,FK] = calsubgray(Index1,ww,SUBindex,np);
%% 3.通过不同方法进行拟合(预计需使用四个不同拟合模型即四个不同函数)
% 1质心法
[ZX1] = ZX(suby,FK);
ZX1 = cell2mat(ZX1);
% 2加权质心法
[ZX2] = JQZX(suby,FK);
ZX2 = cell2mat(ZX2);
% 3高斯拟合
G = 10;
[ZX3] = calGauss(suby,FK,G);
ZX3 = cell2mat(ZX3);
% 4最小二乘拟合
[ZX4] = calOLS(Index1,ww);
ZX4 = cell2mat(ZX4);
%% 图像显示,按照拟合原理,拟合所得点应该位于光斑最亮的像素块内
figure
imshow(uint8(im1));
hold on;
%plot(GCP(:,1),GCP(:,2),'r*');
plot(GCP(:,2),GCP(:,3),'r*');
plot(ZX1(:,1),ZX1(:,2),'b*');
plot(ZX2(:,1),ZX2(:,2),'y*');
plot(ZX3(:,1),ZX3(:,2),'g*');
plot(ZX4(:,1),ZX4(:,2),'m*');
plot(uv(:,1),uv(:,2),'c*');
%%
%%%%存点坐标1
shiftpoint1=  fopen('ZX1.txt','wt'); %fopen('tongmingdian.txt','wt');
for j=1:size(ZX1,1)
    lx=num2str(ZX1(j,1));
    ly=num2str(ZX1(j,2));
    DH = num2str(GCP(j,1));
    fprintf(shiftpoint1,'%c',DH);
    fprintf(shiftpoint1,'%c',' ');
    fprintf(shiftpoint1,'%c',lx);
    fprintf(shiftpoint1,'%c',' ');
    fprintf(shiftpoint1,'%c',ly);
    fprintf(shiftpoint1,'%c\n',' ');
end
fclose(shiftpoint1);
%%
%%%%存点坐标2
shiftpoint2=  fopen('ZX2.txt','wt'); %fopen('tongmingdian.txt','wt');
for j=1:size(ZX2,1)
    lx=num2str(ZX2(j,1));
    ly=num2str(ZX2(j,2));
    DH = num2str(GCP(j,1));
    fprintf(shiftpoint2,'%c',DH);
    fprintf(shiftpoint2,'%c',' ');
    fprintf(shiftpoint2,'%c',lx);
    fprintf(shiftpoint2,'%c',' ');
    fprintf(shiftpoint2,'%c',ly);
    fprintf(shiftpoint2,'%c\n',' ');
end
fclose(shiftpoint2);
% %%
% %%%%存点坐标1
% shiftpoint3=  fopen('ZX3.txt','wt'); %fopen('tongmingdian.txt','wt');
% for j=1:size(ZX3,1)
%     lx=num2str(ZX3(j,1));
%     ly=num2str(ZX3(j,2));
%     DH = num2str(GCP(j,1));
%     fprintf(shiftpoint3,'%c',DH);
%     fprintf(shiftpoint3,'%c',' ');
%     fprintf(shiftpoint3,'%c',lx);
%     fprintf(shiftpoint3,'%c',' ');
%     fprintf(shiftpoint3,'%c',ly);
%     fprintf(shiftpoint3,'%c\n',' ');
% end
% fclose(shiftpoint3);
% %%
% %%%%存点坐标1
% shiftpoint4=  fopen('ZX4.txt','wt'); %fopen('tongmingdian.txt','wt');
% for j=1:size(ZX4,1)
%     lx=num2str(ZX4(j,1));
%     ly=num2str(ZX4(j,2));
%     DH = num2str(GCP(j,1));
%     fprintf(shiftpoint4,'%c',DH);
%     fprintf(shiftpoint4,'%c',' ');
%     fprintf(shiftpoint4,'%c',lx);
%     fprintf(shiftpoint4,'%c',' ');
%     fprintf(shiftpoint4,'%c',ly);
%     fprintf(shiftpoint4,'%c\n',' ');
% end
% fclose(shiftpoint4);
%% 3 精度验证

