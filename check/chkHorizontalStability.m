%% chkStability Each bin ���ĸ�������ȶ����ж�
%% Form
%    chkHorizontalStability
%    ����bin���Ƿ�ˮƽ�ȶ�����Ҫ��������;�е��ȶ��ԣ�����ÿ�����̵������Ƿ����������̻���wall֧��

% ���ڶ��壺ǰ�������ĸ�������LU��maxAdjDist��100���ڵ����������뱾���̵Ľ��������У�����Ϊ�÷������������̣��������������̡�
% NOTE: 
% ?	����ǰ�⣺���������ڣ�
% ?	�Ҳࣺ��������Ϊ˦β����/��������Ҳ���һ����Ϊ˦β���̣���ʹ�����ڡ�
% ?	��ࣺ�����Ƿ����������̣�����Ϊ���ڡ�

% subTLU - ����ĳ��bin�ڵ�����
% subVeh - ĳ������
function isStablility = chkHorizontalStability(subTLU,subVeh)

        global ISplotShowGapAdjust
        isStablility = 1;
        maxAdjDist = 100; 
        
%         subTLU=sortrows(subTLU,{'LU_Bin'},{'ascend'}); 
        subTLU=sortrows(subTLU,{'CoordLUBin'},{'ascend'}); 
        
            subTLU.leftStab =   zeros(height(subTLU),1);
            subTLU.rightStab = zeros(height(subTLU),1);
            subTLU.frontStab = zeros(height(subTLU),1);
            subTLU.backStab =  zeros(height(subTLU),1);
            
        if istable(subTLU),   [subTLU] = getSturctT(subTLU); end
        if istable(subVeh),   [subVeh] = getSturctT(subVeh); end

        % ѭ��ÿ������, �����ж��Ƿ�stability
        for idxl=1:length(subTLU.ID)
            % ����polygon, �����Ŀ�,�� ����pgVEHLeft + �����ĳ�,�� ����pgVEHFront
            pgVEHLeft = polyshape(pgRectangle(0,0,subVeh.LWH(1),subVeh.LWH(3)));
            pgVEHFront = polyshape(pgRectangle(0,0,subVeh.LWH(2),subVeh.LWH(3)));
            
            % ��ǰ���̵�����ֵ
            x=subTLU.CoordLUBin(1,idxl);
            y=subTLU.CoordLUBin(2,idxl);
            z=subTLU.CoordLUBin(3,idxl);
            
            % ��ǰ���̵Ŀ��� - ��Ӧxyz
            w = subTLU.LWH(1,idxl);
            l = subTLU.LWH(2,idxl);
            d = subTLU.LWH(3,idxl);
            
            % ��ǰ���̵�margin (��������)
            m1 = subTLU.margin(1,idxl);
            m2 = subTLU.margin(2,idxl);
            m3 = subTLU.margin(3,idxl);
            m4 = subTLU.margin(4,idxl);
            
            % ����polygon, �����̵Ŀ�,�� ����pgLULeft + �����̵ĳ�,�� ����pgLUFront
            pgLULeft = polyshape(pgRectangle(x,z,w,d));
            pgLUFront = polyshape( pgRectangle(y,z,l,d));
            
            % ����ֵ: �������������̵��������λ�� (���ڵĶ���: ���뱾���̵��������3��margin֮�� �� 30cm֮�� (��ʱmarginΪ0) )
            YleftLU =  ( y - max(3*m1,maxAdjDist));  % ��ǰLU���Yֵ
            XfrontLU =  ( x - max(3*m2,maxAdjDist));  % ��ǰLUǰ��xֵ
            YrightLU =  ( y + l + max(3*m3,maxAdjDist));  % ��ǰLU�Ҳ�Yֵ
            XbackLU =  ( x + w + max(3*m4,maxAdjDist));  % ��ǰLU���xֵ
            
            % ����ֵ: ��������, ��ȷ���Ƿ����ڵ����̼��ϵĶ�Ӧ����ֵ
            Yright = subTLU.CoordLUBin(2,:)+subTLU.LWH(2,:);  % ����LU���Ҳ�Y����ֵ (Yֵ+LU�ĳ�) % ��ǰ������margin��Χ�ڵ��������̣���Ӧ���̵��Ҳࣩ
            Xfront= subTLU.CoordLUBin(1,:);  % ����LU��ǰ��X����ֵ (Xֵ+LU�Ŀ�)
            Yleft = subTLU.CoordLUBin(2,:);  % ����LU�����Y����ֵ (Yֵ)
            Xback = subTLU.CoordLUBin(1,:)+subTLU.LWH(1,:);  % ����LU�ĺ�X����ֵ (Xֵ)
            
            % �ж��ĸ������Ƿ������Ҫ��Horizontal���ȶ���Ҫ��
            [flagLeftStability,flagRightStability,flagFrontStability,flagBackStability]  = deal(zeros(1,1));
            
            %% 1 �ж����
            if YleftLU<=0  % ������ǳ������wall, ���ȶ�
                flagLeftStability = 1;
                pgLUSet =pgVEHLeft;
            else % ����, �ҳ����ص���λ���������Ե����ֵ�ڵ�LU����pgLUSet(����������)
                setidxLU = (Yright >= YleftLU) & (Yright <= y);  setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,1);
                
%                 yc = y - Yright;                yc(yc<0) = Inf;
%                 [a,b] = min(yc);
            end
            % �����̺�������̼��Ľ���. ������, ���ȶ�; ����, ���ȶ�
            pgGap = intersect(pgLUSet,pgLULeft);     if ~isempty(pgGap.Vertices),      flagLeftStability = 1;       end
            
            if flagLeftStability==0 && ISplotShowGapAdjust
                plotStabilityPg(pgVEHLeft,pgLULeft,pgLUSet,pgGap,subVeh.LWH,1);    
                subTLUTTT = struct2table(structfun(@(x) x',subTLU,'UniformOutput',false));
                setidxLUTTT = setidxLU;            setidxLUTTT(idxl) = 1;
                plotSolutionT(subTLUTTT(setidxLUTTT,:),subVeh,0,0,0,1,3,'subTLUUUUUUUU'); 
                plotSolutionT(subTLUTTT,subVeh,0,0,0,1,3,'subTLUUUUUUUU'); 
            end
            

            %% 2 �ж��Ҳ� ��˦β���� �������ȶ���; ���˦β����, �����Ҳ�������˦β���̣�Ҳ�������ȶ���
            if YrightLU >= subVeh.LWH(2) % ���Ҳ��ǳ����ұ�wall, ���ȶ�
                flagRightStability = 1;
                pgLUSet =pgVEHLeft;
            else  % ����, �ҳ����ص���λ������Ҳ��Ե����ֵ�ڵ�LU����pgLUSet(����������)
                setidxLU = (Yleft <= YrightLU) & (Yleft >= y + l);           setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,1);
            end
            % �����̺��Ҳ����̼��Ľ���. ������, ���ȶ�; ����, ���ȶ�
            pgGap = intersect(pgLUSet,pgLULeft);     if ~isempty(pgGap.Vertices),      flagRightStability = 1;       end
            if (flagRightStability==0 && subTLU.isShuaiWei(idxl)==1), flagRightStability=1; end % ���⴦��, right������˦βLU����Ϊunstability.
            if flagRightStability==0 && any(subTLU.isShuaiWei(setidxLU)==1), flagRightStability=1; end % ���⴦��, right������������һ�Ҳ���˦βLU����Ϊunstability.
            

            if flagRightStability==0 && ISplotShowGapAdjust
                plotStabilityPg(pgVEHLeft,pgLULeft,pgLUSet,pgGap,subVeh.LWH,1);           
                subTLUTTT = struct2table(structfun(@(x) x',subTLU,'UniformOutput',false));
                setidxLUTTT = setidxLU;                 setidxLUTTT(idxl) = 1;
                plotSolutionT(subTLUTTT(setidxLUTTT,:),subVeh,0,0,0,1,3,'subTLUUUUUUUU');
                plotSolutionT(subTLUTTT,subVeh,0,0,0,1,3,'subTLUUUUUUUU');
            end
                        
            %% 3 �ж�ǰ��
            if XfrontLU<=0 % ��ǰ���ǳ���ǰ��wall, ���ȶ�
                flagFrontStability = 1;
                pgLUSet =pgVEHFront;
            else  % ����, �ҳ����ص���λ�����ǰ���Ե����ֵ�ڵ�LU����pgLUSet(����������)
                setidxLU = (Xback >= XfrontLU) & (Xback <= x);            setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,2);
            end
            % �����̺�ǰ�����̼��Ľ���. ������, ���ȶ�; ����, ���ȶ�
            pgGap = intersect(pgLUSet,pgLUFront);             if ~isempty(pgGap.Vertices),      flagFrontStability = 1;       end 
            
           if flagFrontStability==0 && ISplotShowGapAdjust
               plotStabilityPg(pgVEHFront,pgLUFront,pgLUSet,pgGap,subVeh.LWH,2);  
               subTLUTTT = struct2table(structfun(@(x) x',subTLU,'UniformOutput',false));
               setidxLUTTT = setidxLU;                setidxLUTTT(idxl) = 1;
               plotSolutionT(subTLUTTT(setidxLUTTT,:),subVeh,0,0,0,1,3,'subTLUUUUUUUU');
               plotSolutionT(subTLUTTT,subVeh,0,0,0,1,3,'subTLUUUUUUUU');
           end
            
                       
            %% 4 �ж���� Back�����Լ��Ķ�Ĭ��Ϊ1,��ʼ�մ���back�ȶ�         
            if XbackLU>= subVeh.LWH(1) % �����ǳ�����wall, ���ȶ�
                flagBackStability = 1;
                pgLUSet =pgVEHFront;
            else  % ����, �ҳ����ص���λ��������Ե����ֵ�ڵ�LU����pgLUSet(����������)
                setidxLU = (Xfront <= XbackLU) & (Xfront >= x + w);                 setidxLU(idxl) = 0;
                pgLUSet = getPgLUSet(setidxLU,subTLU,2);
            end
            % �����̺ͺ�����̼��Ľ���. ������, ���ȶ�; ����, ���ȶ�
            pgGap = intersect(pgLUSet,pgLUFront);             if ~isempty(pgGap.Vertices),      flagBackStability = 1;       end
            if flagBackStability==0, flagBackStability=1; end % ���⴦��, back������Զ����Ϊunstability.
            
            if flagBackStability==0 && ISplotShowGapAdjust % �����ܷ���
                plotStabilityPg(pgVEHFront,pgLUFront,pgLUSet,pgGap,subVeh.LWH,2);                         
            end

            
            %% ��¼ÿ��
            subTLU.leftStab(idxl) =   flagLeftStability;
            subTLU.rightStab(idxl) =   flagRightStability;
            subTLU.frontStab(idxl) =   flagFrontStability;
            subTLU.backStab(idxl) =   flagBackStability;
        end
        
        % �����κ�һ��LU�ǲ�Stability. ��Gap����ʧЧ
        if ~all(subTLU.leftStab) || ~all(subTLU.rightStab) || ~all(subTLU.frontStab) || ~all(subTLU.backStab)
            find(subTLU.leftStab==0);
            find(subTLU.rightStab==0);
            find(subTLU.frontStab==0);
            find(subTLU.backStab==0);
            isStablility = 0;
        end
        
end



% ����װ��϶ͼ����
function plotStabilityPg(pgVEH,pgLU,pgLUSet,pgGap,LWH,type)       
    if type==1
        a = LWH(1);      
    elseif type==2
        a = LWH(2);  
    else
        error('1');
    end
    plot(pgVEH,'FaceColor','white','FaceAlpha',0.01);         hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    plot(pgLU,'FaceColor','green','FaceAlpha',0.2);           hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    plot(pgLUSet,'FaceColor','red','FaceAlpha',0.2);            hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    plot(pgGap,'FaceColor','red','FaceAlpha',0.9);            hold on;    axis equal;    grid on;    xlim([0 1.5*a]);    ylim([0 1.2*LWH(3)]);
    hold off; 

end
  

% ��ȡ���жϵ��ڵ�, ��Boundary��Centroid��ͬ����
function pgLUSet = getPgLUSet(setidxLU,subTLUNewBuff,type);
pgLUSet = polyshape();
if sum(setidxLU) == 0,    return; end
P = [];
setidxLU = find(setidxLU);
for idxl2=1:length(setidxLU)
    fidx2 = setidxLU(idxl2);
    % ���ù̶�����, ���пɱ��������
    x1=subTLUNewBuff.CoordLUBin(1,fidx2);
    y1=subTLUNewBuff.CoordLUBin(2,fidx2);
    z1=subTLUNewBuff.CoordLUBin(3,fidx2);
    
    w1 = subTLUNewBuff.LWH(1,fidx2);
    l1 = subTLUNewBuff.LWH(2,fidx2);
    d1 = subTLUNewBuff.LWH(3,fidx2);
    
    if type == 1
        P =  [P;pgRectangle(x1,z1,w1,d1);[NaN,NaN]];
    elseif type==2
        P =  [P;pgRectangle(y1,z1,l1,d1);[NaN,NaN]];
    else
        error('1');
    end
end
pgLUSet = polyshape(P);
end