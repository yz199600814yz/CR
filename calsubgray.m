function [suby,FKG,FKID,FK] = calsubgray(idgray,gray,subIN,np)
%  此函数通过反距离权重计算亚像素灰度值并进行存储
%  输入参数有idgray,为矩形窗口的索引元胞数组,gray为整数索引处的灰度值,subIN亚像素索引,输入参数都为n*1的元胞数组
%% step1 将元胞数组进行分组分别进行计算
FK = {};
FKG = {};
FKID = {};
for i = 1:size(idgray,1)
    tempid = idgray{i,:};
    tempgray = double(gray{i,:});
    tempsub = subIN{i,:};
    [m,n] = size(tempgray);
    m1 = m-1; n1 = n-1;
    divide = cell(m1,n1);
    %% 对其进行分组
    % step1 首先对灰度整数索引矩阵分组
    tempy = [];
    for i = 1:size(tempid,1)
        tempx = [];
        for j =1:size(tempid,2)
            while j<size(tempid,2)&&i<size(tempid,1)
                temp = cell2mat(tempid(i:i+1,j:j+1));
                tempx = [tempx temp];
                break
            end   
        end
        tempy = [tempy;tempx];
    end
    bb = zeros(1,m1);
    for ii = 1:size(bb,2)
        bb(ii) = 1;
    end
    dim1 = bb.*2;
    dim2 = bb.*4;
    tempy = mat2cell(tempy,dim1,dim2);
    
    % step2 对矩形窗口灰度进行分组
    tempG = [];
    for i = 1:size(tempgray,1)
        tempGX = [];
        for j =1:size(tempgray,2)
            while j<size(tempgray,2)&&i<size(tempgray,1)
                temp = tempgray(i:i+1,j:j+1);
                tempGX = [tempGX temp];
                break
            end   
        end
        tempG = [tempG;tempGX];
    end
    dim3 = bb.*2;
    dim4 = bb.*2;
    tempG =  mat2cell(tempG,dim3,dim4);
    
    
    % step3 对亚像素索引进行分组
    % (1) 首先找到uv均为整数像素行列号
    % (2) 使用整数索引对其进行分割   
    suby = {};
    for i = 1:np:size(tempsub,1)
        subx = [];
        for j = 1:np:size(tempsub,2)
            while j<size(tempsub,2)&&i<size(tempsub,1)
                temp = {tempsub(i:i+np,j:j+np)};
                subx = [subx temp];
                break
            end
        end
        suby = [suby;subx];
    end
    tempm1 = {tempy};
    tempm2 = {tempG};
    tempm3 = {suby};
    FKG = [FKG;tempm2];      % 分块矩阵灰度
    FKID = [FKID;tempm1];    % 分块矩阵索引
    FK = [FK;tempm3];             % 亚像素分块索引
end
%% step2 对分块后的亚像素进行灰度计算
suby = {};
for i = 1:size(idgray,1)
    tempfkg = FKG(i,:);
    tempsub = FK(i,:);
    tempfkid = FKID(i,:);
    for j = 1:size(tempfkg,1)
        for k = 1:size(tempfkg,2)
            aa = tempfkg{j,k};
            bb = tempsub{j,k};
            cc = tempfkid{j,k};
            yyyyyy = {};
            for ii = 1:size(aa,1)
                xxxxxx = {};
                for jj = 1:size(aa,2)
                    aaa = aa{ii,jj};
                    bbb = bb{ii,jj};
                    bbb(1,:) = [];
                    bbb(:,1) = [];
                    ccc = cc{ii,jj};
                    tempbbby = [];
                    for iii = 1:size(bbb,1)
                        tempbbbx =[];
                        for jjj = 1:size(bbb,2)
                            tempbbb = bbb{iii,jjj};
                            subgrayzs = ((1/sqrt((tempbbb(:,1)- ccc(1,1))^2+(tempbbb(:,2)-ccc(1,2))^2))^0.5)*aaa(1,1);
                            subgrayys = ((1/sqrt((tempbbb(:,1)- ccc(1,3))^2+(tempbbb(:,2)-ccc(1,4) )^2))^0.5)*aaa(1,2);
                            subgrayzx = ((1/sqrt((tempbbb(:,1)- ccc(2,1))^2+(tempbbb(:,2)-ccc(2,2) )^2))^0.5)*aaa(2,1);
                            subgrayyx = ((1/sqrt((tempbbb(:,1)- ccc(2,3))^2+(tempbbb(:,2)-ccc(2,4) )^2))^0.5)*aaa(2,2);                                        
                            subgray = subgrayzs+ subgrayys+subgrayzx+subgrayyx;
                            tempbbbx = [ tempbbbx subgray];
                        end
                        tempbbby = [tempbbby;tempbbbx]; 
                    end
                    tempbbby(np,np) = (tempbbby(np,np-1) +tempbbby(np-1,np))/2 ;
                    xxxxxx = [ xxxxxx tempbbby];
                end
                yyyyyy = [yyyyyy;xxxxxx];
            end
        end
    end
    suby = [suby;{yyyyyy}];
end
end

