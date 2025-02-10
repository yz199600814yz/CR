function [ZXECNHD] = calOLS(suby,FK)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
%% 数据传入
% suby = load('suby.mat');
% suby = suby.suby;   % 分块灰度
% FK = load('FK.mat');
% FK = FK.FK;        % 分块索引
%% 数据合并
zhx = {};
zhy = {};
for i = 1:size(suby,1)
    Ctemp = suby{i,:};     % Ctemp就是cell(temp)
    tempFK = FK{i,:};
    xY = [];
    yY = [];
    for j = 1:size(Ctemp,1)
        xX = [];
        yX = [];
        for k = 1:size(Ctemp,2)
            temp = Ctemp{j,k};
            x0 = temp(:,1);
            y0 = temp(:,2);
            xX = [xX x0];
            yX = [yX y0];
        end
        xY = [xY;xX];
        yY = [yY;yX];
    end
    zhx = [zhx;xY];
    zhy = [zhy;yY];
end
%% 最小二乘计算曲面峰值点
ZXECNHD = {};
for i = 1:size(FK,1)
    tempzhsuby = FK{i,:};
    tempzhx = zhx{i,:};
    tempzhy = zhy{i,:};
    % step1 首先转换为n*1的数据形式
    numP = numel(tempzhsuby);
    tempzhsuby  = reshape(tempzhsuby,numP,1);
    tempzhx = reshape(tempzhx,numP,1);
    tempzhy = reshape(tempzhy,numP,1);
    
    
    A = zeros(6);
    C = zeros(6,1);
    
    for j = 1:size(tempzhx,1)
        tempA = zeros(6);
        tempC = zeros(6,1);
        tempA(1,1) = (tempzhx(j))^4;
        tempA(1,2) = ((tempzhx(j))^2)*((tempzhy(j))^2);
        tempA(1,3) = ((tempzhx(j))^3)*tempzhy(j);
        tempA(1,4) = ((tempzhx(j))^3);
        tempA(1,5) = ((tempzhx(j))^2)*tempzhy(j);
        tempA(1,6) = ((tempzhx(j))^2);
        
        tempA(2,1) = tempA(1,2);
        tempA(2,2) = (tempzhy(j))^4;
        tempA(2,3) = tempzhx(j)*((tempzhy(j))^3);
        tempA(2,4) = tempzhx(j)*((tempzhy(j))^2);
        tempA(2,5) = ((tempzhy(j))^3);
        tempA(2,6) = ((tempzhy(j))^2);
        
        tempA(3,1) = tempA(1,3);
        tempA(3,2) = tempA(2,3);
        tempA(3,3) = tempA(1,2);
        tempA(3,4) = tempA(1,5);
        tempA(3,5) = tempA(2,4);
        tempA(3,6) = tempzhx(j)*tempzhy(j);
        
        tempA(4,1) = tempA(1,4);
        tempA(4,2) = tempA(2,4);
        tempA(4,3) = tempA(3,4);
        tempA(4,4) = tempA(1,6);
        tempA(4,5) = tempA(3,6);
        tempA(4,6) = tempzhx(j);
        
        tempA(5,1) = tempA(1,5);
        tempA(5,2) = tempA(2,5);
        tempA(5,3) = tempA(3,5);
        tempA(5,4) = tempA(4,5);
        tempA(5,5) = tempA(2,6);
        tempA(5,6) = tempzhy(j);
        
        tempA(6,1) = tempA(1,6);
        tempA(6,2) = tempA(2,6);
        tempA(6,3) = tempA(3,6);
        tempA(6,4) = tempA(4,6);
        tempA(6,5) = tempA(5,6);
        A = A+tempA;
        
        tempC(1,:) = ((tempzhx(j))^2)*tempzhsuby(j);
        tempC(2,:) = ((tempzhy(j))^2)*tempzhsuby(j);
        tempC(3,:) = tempzhx(j)*tempzhy(j)*tempzhsuby(j);
        tempC(4,:) = tempzhx(j)*tempzhsuby(j);
        tempC(5,:) = tempzhy(j)*tempzhsuby(j);
        tempC(6,:) = tempzhsuby(j);
        C = C+tempC;
        
    end
    A(6,6) = numP;
    B = A\C;
    xx = ((2*B(2)*B(4)) - (B(3)*B(5)) )/( (B(3)^2) - 4*B(1)*B(2));
    yy = ((2*B(1)*B(5)) - (B(3)*B(4)) )/( (B(3)^2) - 4*B(1)*B(2));
    tempNH = [xx yy];
    ZXECNHD = [ZXECNHD;tempNH];
end
end

