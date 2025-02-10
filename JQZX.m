function [ZX] = JQZX(suby,FK)
% ��Ȩ�����㷨
%   �˴���ʾ��ϸ˵��
%  ͨ�����ķ��������
%   ͨ������ÿһ��ϸ�����ص��������Ӧ�Ҷ�ֵ�ĳ˻������Ӧ���ؿ������
%% ���ݴ���
% suby = load('suby.mat');
% suby = suby.suby;   % �ֿ�Ҷ�
% FK = load('FK.mat');
% FK = FK.FK;        % �ֿ�����
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
%% ��Ȩ���ķ�����
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

