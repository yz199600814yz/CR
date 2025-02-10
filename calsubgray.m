function [suby,FKG,FKID,FK] = calsubgray(idgray,gray,subIN,np)
%  �˺���ͨ��������Ȩ�ؼ��������ػҶ�ֵ�����д洢
%  ���������idgray,Ϊ���δ��ڵ�����Ԫ������,grayΪ�����������ĻҶ�ֵ,subIN����������,���������Ϊn*1��Ԫ������
%% step1 ��Ԫ��������з���ֱ���м���
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
    %% ������з���
    % step1 ���ȶԻҶ����������������
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
    
    % step2 �Ծ��δ��ڻҶȽ��з���
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
    
    
    % step3 ���������������з���
    % (1) �����ҵ�uv��Ϊ�����������к�
    % (2) ʹ����������������зָ�   
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
    FKG = [FKG;tempm2];      % �ֿ����Ҷ�
    FKID = [FKID;tempm1];    % �ֿ��������
    FK = [FK;tempm3];             % �����طֿ�����
end
%% step2 �Էֿ��������ؽ��лҶȼ���
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

