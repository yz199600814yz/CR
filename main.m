clc;
clear;
close all;
%% ����׼����Ӱ��Ϳ��Ƶ��ȡ
% im1 = uint8(imread('Ӱ��\ZW-3M-5624.tiff'));
im1 = imread('Ӱ��\SS-1m-2170.tiff');
GCP = load('���Ƶ�\SS-1m-2170.txt');
% GCP = GCP(:,2:3);
pointnum = size(GCP,1);
%% ��ȡ���Ƶ���δ���
% step1 ��������ת��Ϊ�������겢ȡ�����ڵĻҶ�ֵ����������
win = 1;               % ���ڴ�С��ֵ
ww = cell(pointnum,1); % �������������ڴ洢���������ĻҶ�ֵ
Index1 = cell(pointnum,1); %�������������ڴ洢��������ֵ
uv = [];                   %���Ƶ���������
for i = 1:pointnum 
    %temp = make_u_v(GCP(i,1),GCP(i,2));
    temp = make_u_v(GCP(i,2),GCP(i,3));
    tempw = im1(temp(:,2)-win:temp(:,2)+win,temp(:,1)-win:temp(:,1)+win);
    % ȡ��ȡ���Ƶ�ľ�����������
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
    ww{i,:} = tempw;   % ���ο�Ҷ�Ԫ������
    Index1{i,:} = Index;  %���ο���������Ԫ������
    uv = [uv;temp];
end
%%  �Դ��������ؽ��������طָ�
% step2 ȡ�������ĵ㲢�趨���������ؽ��зָͨ��Ȩ�غ����Էָ����ؿ���и�ֵ
np = 10;              %�趨�ָ���
length = 1/np;                 %�ָ��
% �Ծ��󴰿��ڵ����ؽ��������طָ�
% 1.ȡԪ����������ϽǺ����½���
SUBindex = cell(pointnum,1);
for i = 1:size(Index1,1)
    winpt = Index1{i,:};             %3*3��������Ԫ������
    [m,n] = size(winpt);
    minuv = winpt{1,1};
    maxuv = winpt{m,n};
    x0 = [minuv(:,1):length:maxuv(:,1)];
    y0 = [minuv(:,2):length:maxuv(:,2)];
    tempX = [];
    for j = 1:size(x0,2)
        tempY = [];
        for k = 1:size(y0,2)
            % ����������ϸ��
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
%% 2. ͨ��Ȩ�غ�������ÿ����ĻҶ�ֵ
[suby,FKG,FKID,FK] = calsubgray(Index1,ww,SUBindex,np);
%% 3.ͨ����ͬ�����������(Ԥ����ʹ���ĸ���ͬ���ģ�ͼ��ĸ���ͬ����)
% 1���ķ�
[ZX1] = ZX(suby,FK);
ZX1 = cell2mat(ZX1);
% 2��Ȩ���ķ�
[ZX2] = JQZX(suby,FK);
ZX2 = cell2mat(ZX2);
% 3��˹���
G = 10;
[ZX3] = calGauss(suby,FK,G);
ZX3 = cell2mat(ZX3);
% 4��С�������
[ZX4] = calOLS(Index1,ww);
ZX4 = cell2mat(ZX4);
%% ͼ����ʾ,�������ԭ��,������õ�Ӧ��λ�ڹ�����������ؿ���
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
%%%%�������1
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
%%%%�������2
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
% %%%%�������1
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
% %%%%�������1
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
%% 3 ������֤

