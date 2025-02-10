function [GSNHD] = calGauss(suby,FK,G)
%   ʹ�ö�ά��˹�����������������
%   �˴���ʾ��ϸ˵��
% ���ݴ���
% suby = load('suby.mat');
% suby = suby.suby;   % �ֿ�Ҷ�
% FK = load('FK.mat');
% FK = FK.FK;        % �ֿ�����
% G = 10;
%% ���ݺϲ�
zhsuby = {};
zhx = {};
zhy = {};
for i = 1:size(suby,1)
    Ctemp = suby{i,:};     % Ctemp����cell(temp)
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
%%  ��˹���
% syms x0 y0 a b    % ����δ֪�� x0,y0Ϊ����ϵ�����,aΪx�����׼��,bΪy�����׼��
GSNHD = {};
for i = 1:size(zhsuby,1)
    tempzhsuby = zhsuby{i,:};
    tempzhx = zhx{i,:};
    tempzhy = zhy{i,:};
    % step1 ����ת��Ϊn*1��������ʽ
    numP = numel(tempzhsuby);
    tempzhsuby  = reshape(tempzhsuby,numP,1);
    tempzhx = reshape(tempzhx,numP,1);
    tempzhy = reshape(tempzhy,numP,1);
    A = [];
    B = [];
    for j = 1:size(tempzhx,1)
%         c=(tempzhx(j) - x0)^2;
%         d=(tempzhy(j) - y0)^2;
%         e=1/(2*a^2); 
%         f= 1/(2*b^2);
%         lnF = (log(G) - (x0^2)*e - (y0^2)*f )+x0*2*e*tempzhx(j) + y0*2*f*tempzhy(j) - e*(tempzhx(j))^2 - f*(tempzhy(j))^2;
        ai = tempzhsuby(j)*log(tempzhsuby(j));
        bi = [tempzhsuby(j) tempzhsuby(j)*tempzhx(j) tempzhsuby(j)*tempzhy(j) tempzhsuby(j)*(tempzhx(j))^2 tempzhsuby(j)*(tempzhy(j))^2];
        A = [A;ai];
        B = [B;bi];
    end
%     C = [(log(G) - (x0^2)*e - (y0^2)*f ) (x0*2*e) (y0*2*f) (-1*e) (-1*f)]';
    [Q,R] = qr(B);
    AA = Q'*A;
    S = AA(1:5,:);
    R1 = R(1:5,:);
    CC = R1\S;
    xx = (-1/2)*(CC(2)/CC(4));
    yy = (-1/2)*(CC(3)/CC(5));
    tempNH = [xx yy];
    GSNHD = [GSNHD;tempNH];
end 
end

