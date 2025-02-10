function [ZX] = JQZX(suby,FK)
% 加权质心算法
%   此处显示详细说明
%  通过质心法进行拟合
%   通过计算每一个细分像素的索引与对应灰度值的乘积计算对应像素块的质心
%% 数据传入
% suby = load('suby.mat');
% suby = suby.suby;   % 分块灰度
% FK = load('FK.mat');
% FK = FK.FK;        % 分块索引
%% 数据合并
zhsuby = {};
zhx = {};
zhy = {};
for i = 1:size(suby,1)
    Ctemp = suby{i,:};     % Ctemp就是cell(temp)
    Mtemp = cell2mat(Ctemp);
    tempFK = FK{i,:};
    xxxxxY = [];
    yyyyyY = [];
    for j = 1:size(Ctemp,1)
        xxxxxX = [];
        yyyyyX = [];
        for k = 1:size(Ctemp,2)
            CFK = tempFK{j,k};
            CFK(1,:) = [];
            CFK(:,1) = [];
            Xy = [];
            Yy = [];
            for jj = 1:size(CFK,1)
                 Xx = [];
                 Yx =[];
                for kk = 1:size(CFK,2)
                    x0 = CFK{jj,kk}(:,1);
                    y0 = CFK{jj,kk}(:,2);
                    Xx = [Xx x0];
                    Yx = [Yx y0];
                end
                Xy = [Xy;Xx];
                Yy = [Yy;Yx];
            end
            xxxxxX = [xxxxxX  Xy];
            yyyyyX = [yyyyyX  Yy];
        end 
        xxxxxY = [xxxxxY;xxxxxX];
        yyyyyY = [yyyyyY;yyyyyX];
    end
    zhsuby = [zhsuby;Mtemp];
    zhx = [zhx;xxxxxY];
    zhy = [zhy;yyyyyY];
end
%% 加权质心法计算
ZX = {};
for i = 1:size(zhsuby,1)
    aa = zhsuby{i,:};
    bb = zhx {i,:};
    cc = zhy{i,:};
    xxY = [];
    yyY = [];
    for j = 1:size(aa,1)
        xxX = [];
        yyX = [];
        for k = 1:size(aa,2)
            tempx = (aa(j,k))^2*bb(j,k);
            tempy = (aa(j,k))^2*cc(j,k);
            xxX = [xxX tempx];
            yyX = [yyX tempy];
        end
        xxY = [xxY;xxX];
        yyY = [yyY;yyX];
    end
    aaa = aa.*aa;
    tempzxX = sum(xxY(:))/sum(aaa(:));
    tempzxY = sum(yyY(:))/sum(aaa(:));
    tempm = [tempzxX tempzxY];
    ZX = [ZX;tempm];
end
end

