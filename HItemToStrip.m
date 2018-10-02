% function [d] = HItemToStrip(d,p)
function [LU,Item,Strip]= HItemToStrip(LU,Item,Veh,p)
% ��Ҫ����:Item����Strip�� %  ����:�����(row);  ����:��������(coloum);
% Input ---  ITEM:  ID LWH Weight
% Output --- ITEM: itemorder Item_Strip itemRotaFlag CoordItemStrip
% Output --- Strip: LW Weight
% Item (1 LWH (��֪)
% Item (2 Item_Strip  (dim1:���item��ĳ��strip dim2:item����˳��(��->��) 
% Item (3 CoordItemStrip Item��strip������) 
% Item (4 itemo rder ������Item����˳��)
% Strip (1 LW )
% tmpStrip_Item         (2,n): % ��1��ÿ��Strip�ڵ�Item���� �� ��2��ÿ��Strip�ڵĲ�ͬLUID����
%% ��ʼ��
% nDim Itemά��(2) nItem Item���� nStrip Strip���� 
% widthStrip Strip�����
nDim = size(Item.LWH,1); if nDim ==3, nDim = nDim-1;end
sz = size(Item.LWH);
nItem = sz(2);
wStrip = Veh.LWH(1,1);

%% ���ж�ITEM����Horizontal/Vertical ��ʽ�ڷţ������Ƿ���ת�������жϽ����㷨��˳��
% �����Ƿ�������ת, ֻ���Ƿ���Ҫ��Horizontal/Vertical��ʽ�ڷ�

%  [ItemLWRota, ItemRotaed] = placeItemHori(Item,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% ���ԭ�ⲻ���ķ���ֵ:����ʼֵ
% x = Item.Rotaed
% [Item.Rotaed] = placeItemHori(Item.LWH,Item.isRota,2);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
%         if any(x~=Item.Rotaed),             error('111111111111');         end
% Item.LWH = getRotaedLWH(Item.LWH, Item.Rotaed, LU.buff); 

    %% ITEM���� 555
    % ��ȡItem��˳�� % ITEM��������ʽ �߶�/��̱�
%     printstruct(Item)    
    [Item.itemorder] = getITEMorder(Item,p.whichSortItemOrder );
    
    % ��ȡ��order������ITEM: sItem
    if isSameCol(Item)
        sItem = structfun(@(x) x(:,Item.itemorder),Item,'UniformOutput',false);
    else
        error('����ʹ��structfun');
    end
%     printstruct(sItem)

    % ������;������Ҫԭ��������Gpreproc���Լ�����H/V���ô�����
% %     % 1��2���ڵ�, sortedItemArray��Ӧ��LWHRota��Rotaed������->���󷵻ص�ԭ����ItemArry��
% %     if p.whichRotationHori == 1 % �����ĸ�level,������horizontally��ʽ�ڷ�
% %          x = sItem.Rotaed;
% %          sItem.LWH
% %         [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,1);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��       
% % %          if any(x~=sItem.Rotaed),                 error('111111111111');         end    
% %     end
% %     if p.whichRotationHori == 2 % �����ĸ�level,������vertical��ʽ�ڷ�
% %                     x = sItem.Rotaed
% %         [ sItem.Rotaed] = placeItemHori(sItem.LWH,sItem.isRota,0);  %�ڶ���������1: Hori; 0: Vert������: ԭ�ⲻ��
% % %                     if any(x~=sItem.Rotaed),                                  error('111111111111');         end
% %     end
% %     sItem.LWH = getRotaedLWH(sItem.LWH, sItem.Rotaed, LU.buff); 
% %     sItem.LWH
%% 55 LU->Item->Stripת�� 
% ��ȡ���� �˴�ֻʹ��LWHRota��Rotaed; ��ʹ��LWH
% ItemLWRotaSort = sItem.LWH(1:2,:); %ItemLWSortHori
% ItemRotaedSort = sItem.Rotaed; % ItemRotaSortHori % itemRotaSort = zeros(1,size(Item.LWH,2));
% ItemisRotaSort = sItem.isRota;
% ItemWeightSort = sItem.Weight;

% Itemorder = Item.itemorder;

%%
% 1 Strip��ʼ��
% ��ȡStrip.LW:  �����ɵ�Strip�ĳ����
Strip.LW = zeros(2,nItem);   %strip���� dim2-����(����ߵļ���) (�߶Ƚ����ο�,current �߶�)
Strip.LW(1,:) = wStrip;   %dim1-���ʣ�� 
Strip.Weight = zeros(1,nItem); % ��ʼ��ֵ
    
% 2  ��ʱ
tmpStrip_Item = zeros(2,nItem);  % ��1��ÿ��Strip�ڵ�Item���� �� ��2��ÿ��Strip�ڵĲ�ͬLUID����
% 3 sItem����
% ��ȡCoordItemStripSort  Item��strip������ֵ
% ��ȡItem_Strip: ÿ�������Item���ĸ�Strip��  �Լ�˳��
sItem.Item_Strip = zeros(2,nItem); %dim1:���ڵڼ���level dim2:���ڸ�level�ڼ����ŷ� 555
sItem.CoordItemStrip = zeros(2,nItem); %Item��strip������ֵ

% 55 ��ȡthisLevel - ��ǰitemҪ�����level���
% ѭ����strip�а���item,���̶�item,�仯ѡ��ͬlevel(thisLevel)
% ע�ͣ���ȡ FLAG        �ɷ��µ�ǰitem(iItem)������һ��level�ļ��ϣ� 
% ע�ͣ���ȡ thisLevel   ��FLAG���ҵ���������Ǹ�thisLevel, ��ִ�� insert����

iLevel = 1; iItem = 1; %iStrip����itemʵ��

while 1
    if iItem > nItem, break; end

    % ���ݲ�ͬ�����ҵ��ܷ��뵱ǰitem��strips/levels�е�һ��    
    [thisLevel,iLevel] = getThisLevel(iItem,iLevel,sItem, Strip, p);     %iLevel���ڴκ����ڲ��ϵ�������Զָʾ��ǰ���µ�level
    
    insertItemToStrip(thisLevel,iItem);
    %     plot2DStrip(); %������ͼ    
    iItem = iItem + 1;
end
    
% figure(randi(1000));
% plot2DStrip();  %����������: һ���Ի�ͼ
figure(randi(1000));
plot3DStrip();

% Item��أ����µİ�˳�򷵻أ��޸��µĲ��践�أ�
% LU�ڲ�����,sLU����order�仯����(��ҪΪ��sLU�������ļ�������,Ҫ��˳��ת����)
% ��ȡRotaed : ÿ��item�Ƿ�Rotation�ı�־
% ��ȡLWHRota��ÿ��item���Rotaed��־���õ�LWH
if isSameCol(sItem)
    Item = reorderStruct(Item.itemorder, sItem);
else
    error('����ʹ��structfun');
end

                    % ��ȡItem_Strip : ÿ��Item���ĸ�Strip��  �Լ�˳��
                    % ��ȡCoordItemStrip : ÿ��Item��Strip������
                    %     Item.Item_Strip(:,Item.itemorder) = sItem.Item_Strip;
                    %     Item.CoordItemStrip(:,Item.itemorder) = sItem.CoordItemStrip;    
                    % ItemArray��ת���
                    %         Item.Rotaed(:,Item.itemorder) = sItem.Rotaed;
                    %     ItemLWRota(:,Item.itemorder) = sItem.LWH(1:2,:);
                    %     Item.LWH = [ItemLWRota; Item.LWH(3,:)];             % ����ԭʼ˳�����ת���ItemArray


                    % LUArray��ת���,��ʱ����    ��Ƕ�뵽rotateItem
                %     nbItem=length(Item.Rotaed);
                %     % ѭ��ÿ��item
                %     for idxItem=1:nbItem
                %         tmpflagThisItem = (LU.LU_Item(1,:)==idxItem );
                %         % ��Ӧλ��LU.Rotaed����
                %         if Item.Rotaed(idxItem)
                %             LU.Rotaed(tmpflagThisItem) = ~LU.Rotaed(tmpflagThisItem);
                %             % ��Ӧλ��LU.LWH����
                %             LU.LWH(1, tmpflagThisItem) = Item.LWH(1, idxItem);
                %             LU.LWH(2, tmpflagThisItem) = Item.LWH(2, idxItem);
                %         end
                %     end
    
% Strip���: ��˳�����(ȥ����ʼ���Ķ�����)
% Strip.LW
% ��ȡLWStrip:  �����ɵ�strip�ĳ���
% ��ȡStripWeight:  �����ɵ�strip������
% ���Strip������ȫ����ͬ
if isSameCol(Strip)
    Strip = structfun(@(x) x( : , Strip.Weight(1,:)>0 ), Strip, 'UniformOutput', false);
else
    error('����ʹ��structfun');
end
                 

%% ����script
    % �����Ҫ���:���ÿ��level������ 
%     printscript();
%     printstruct(d);
    
    %% Ƕ�׺���        
    function insertItemToStrip(thisLevel,iItem)
%         %Ϊ��matlab��coder ����
%         sItem.CoordItemStrip=sItem.CoordItemStrip;Strip.LW=Strip.LW;
%         ItemRotaSort=ItemRotaSort;ItemLWSort=ItemLWSort;
%         Strip_Item=Strip_Item;sItem.Item_Strip=sItem.Item_Strip;
        
        % 1 ����Item���Sort����
        %  1.1 ����CoordItemStripSort
        sItem.CoordItemStrip(1,iItem) = wStrip - Strip.LW(1,thisLevel);        %����x����
        sItem.CoordItemStrip(2,iItem) = sum(Strip.LW(2,1:thisLevel-1));      %����y���� %���iLevel=1,�����ߣ�����Ϊ0������Ϊ���
        
        % 2 ����Strip������ݣ�������
        %  2.1 ����LWStrip
        updateLWStrip(); %�����ж��Ƿ�������ת�������ж��Ƿ�������Level�������ж��Ƿ�ǰ�ɷ���
                % % %         if sItem.isRota(iItem) == 1 %��Item������ת
                % % %             % �ж����: 2��
                % % %             isflagCurr = Strip.LW(1,thisLevel) >=  sItem.LWH(1,iItem); %�ж��Ƿ�current's stripʣ���� >= ��ǰ�߶ȣ�����ת��
                % % %             isNewLevel = Strip.LW(1,thisLevel) == wStrip; % �ж��Ƿ� new Level            
                % % %             % ����strip��Ϣ
                % % %             if isNewLevel %�������,�����Է���,���ۺ��ְڷ�,���:ֱ�Ӹ��£�Item�Ѱ�Hori/Vert�ڷŹ���
                % % %                     updateLWStrip();
                % % %             else % �������level ����԰ڷ�,����; ����,����������ת�������
                % % %                     if isflagCurr
                % % %                         updateLWStrip();
                % % %                     else
                % % %                         error('11111111');
                % % % %                         rotateItem(); %����ֱ��ע��
                % % % %                         updateLWStrip();
                % % %                     end
                % % %              end            
                % % %         elseif sItem.isRota(iItem) == 0 %��Item��������ת
                % % %             % ����strip��Ϣ
                % % %             updateLWStrip();
                % % %         end
        
        %  2.2 ����Strip.Strip_Item ��1 ��strip�ڰ�������Item
        tmpStrip_Item(1,thisLevel) = tmpStrip_Item(1,thisLevel) + 1; %ֻҪ��level����һ��item,����������1
        
        %  2.3 ���±�level��Ӧ��StripWeight: 
        Strip.Weight(thisLevel) =  Strip.Weight(thisLevel) + sItem.Weight(iItem);
        
        %  1.3 ����item����strip��ϢitemBeStripMatrixSort
        sItem.Item_Strip(1,iItem) = thisLevel;    %�ڼ���level
        sItem.Item_Strip(2,iItem) = tmpStrip_Item(1,thisLevel); %��level�µڼ��ΰ���
        
                        %  2.2 ����Strip.Strip_Item ��2 ��strip�ڰ�������Item
                %         itemThisLevel = sItem.Item_Strip(1,:) == thisLevel;
                %         Strip.Strip_Item(2,thisLevel) = numel(unique(sItem.LID(1,itemThisLevel)));

        % 4 ����Ƕ�׺���- ����������Ҫ
% %         function rotateItem()
% %             %  1 �������Rotaed�仯 ��Ҫ��ITEM������rotate(��)��ȥ
% %             sItem.Rotaed(iItem) = ~sItem.Rotaed(iItem);
% %             tep = sItem.LWH(1,iItem);
% %             sItem.LWH(1,iItem) = sItem.LWH(2,iItem);
% %             sItem.LWH(2,iItem) = tep;
% %             sItem.LWH(3,iItem) = sItem.LWH(3,iItem);
% %             
% %             %  2 Item�ڲ���LuҲҪ��֮����
% %             tmpflagThisItem = (LU.LU_Item(1,:)==iItem);
% %             % ��Ӧλ��LU.Rotaed����
% %             if Item.Rotaed(iItem)
% %                 LU.Rotaed(tmpflagThisItem) = ~LU.Rotaed(tmpflagThisItem);
% %                 % ��Ӧλ��LU.LWH�ĳ�����£��߶Ȳ�����
% %                 LU.LWH(1, tmpflagThisItem) = Item.LWH(1, iItem);
% %                 LU.LWH(2, tmpflagThisItem) = Item.LWH(2, iItem);
% %             end
% %         end
        
        function updateLWStrip()
            Strip.LW(1,thisLevel) = Strip.LW(1,thisLevel) - sItem.LWH(1,iItem); %����wleft (�ڷŷ���ǰ��һ��)
            Strip.LW(2,thisLevel) = max(Strip.LW(2,thisLevel), sItem.LWH(2,iItem)); %����strip�߶�lleft(ȡ���ֵ)
            
            % ����Strip�а���ID�����
%             Strip.LID(sItem.LID(1,iItem),thisLevel) =  1;         
%             Strip.SID(sItem.SID(1,iItem),thisLevel) =  1;        
%             Strip.UID(sItem.UID(1,iItem),thisLevel) = 1;         % ��ֵΪ�������

%              Strip.PID(:,thisLevel) = Strip.PID(:,thisLevel) + sItem.PID(:,iItem); % ��ֵΪ���ִ���
%              Strip.PID(Strip.PID>0) = 1; % ��ֵ��Ϊ�������
            
%             Strip.SID(sItem.SID(1,iItem),thisLevel) = 1;         % 555 ���¶���PID
%             Strip.UID(sItem.UID(1,iItem),thisLevel) = 1;         % 555 ���¶���PID
%             Item.PID(sLU.PID(1,iLU),thisItem) = 1;         % 555 ���¶���PID
            
        end
        
    end
    
    function printscript()
        % ���Դ���
        % % LWHStrip
        % % sItem.Item_Strip
        % % ItemLWSort
        % % itemCoordMatrixSort
        % % Item_Strip
        % % LWHItem
        % % itemCoordMatrix
        %  printstruct(d);
        
        % �����Ҫ���:��ô�1��ʼÿ��strip����������
        for iStrip = 1:max(Item.Item_Strip(1,:))
            [~,idx] = find(Item.Item_Strip(1,:)==iStrip);
            fprintf('strip %d ��ʣ���+���Ϊ:  ',iStrip);
            fprintf('( %d ) ',Strip.LW(:,iStrip));
            fprintf('\n');
            fprintf('strip %d ���� original Item ������(����)[��ת��־]{����}Ϊ  \n  ',iStrip);
            fprintf('%d ',idx);
            fprintf('( %d ) ', Item.LWH(1:nDim,idx));fprintf('\n');
            fprintf('[ %d ] ', Item.Rotaed(:,idx));fprintf('\n');  %sItem.Rotaed
            fprintf('{ %d } ', Item.CoordItemStrip(:,idx));fprintf('\n');
            fprintf('\n');
        end
    end
    
    % sItem������ά��ͼ
    function plot3DStrip()
        nIDType = unique(LU.ID);
        nColors = hsv(length(nIDType)); %��ͬ����Item��LU���費ͬ��ɫ
        
        yxz = sItem.LWH;
        coord = sItem.CoordItemStrip;   coord(3,:) = 0;
        
        for i=1:numel(sItem.Weight)
            j=sItem.LID(i);
            LUColor = 0.8*nColors(nIDType==j{1}, : );
            plotcube(yxz(:,i)',coord(:,i)',0.7,LUColor);
            
            % Set the lable and the font size
            axis equal;         grid on;        view(60,40);
            xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
            xlim([0 Veh.LWH(1,1)]);         zlim([0 Veh.LWH(3,1)]);
        end
    end

    function plot2DStrip()
        %% ��ʼ��
        wStrip = wStrip;        
        hStrip = sum(Strip.LW(2,sItem.Item_Strip(2,:)>0));        
        nstrip = sum(sItem.Item_Strip(2,:)>0);

        tmpLID = cellfun(@(x)x(1), sItem.LID);        
        nIDType = unique(tmpLID);
        nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
        
        %% ��ͼ
        % 1 ��ͼ��������Strip
        DrawRectangle([wStrip/2 hStrip/2 wStrip hStrip 0],'--', [0.5 0.5 0.5]);
        hold on;
        % 2 ��ͼ�����strip/item ��ͼ
        for istrip = 1:nstrip
            % �ҳ���ǰistrip����Ʒ����
            idxDrawItem = find(sItem.Item_Strip(1,:)==istrip);
            % ��ȡ�������µı���
            drawItemCoordMatrix = sItem.CoordItemStrip(:,idxDrawItem);
            drawItemLWH = sItem.LWH(:,idxDrawItem);
            drawItemId = tmpLID(:,idxDrawItem);

            % ��ͼ�����item
            nThisItem = size(drawItemLWH,2);
            for iplotItem = 1:nThisItem
                % ��ͼ��������iItem
                itemWidth = drawItemLWH(1,iplotItem);
                itemLength = drawItemLWH(2,iplotItem);
                itemCenter = [drawItemCoordMatrix(1,iplotItem)+itemWidth/2 ...
                    drawItemCoordMatrix(2,iplotItem)+itemLength/2 ];
                
                % ���ӶԱ���iItem�����ͣ���ɫ���ж�
                itemID = drawItemId(iplotItem);
                itemColor = 0.8*nColors(nIDType==itemID, : );
                
                DrawRectangle([itemCenter itemWidth itemLength 0],  '-',itemColor);
                hold on;
            end
        end
        % hold off;
    end
end

%% �ֲ����� %%

%% ����1: getITEMorder
% ����ITEM��˳��,��NEXT FIT�ķ�ʽ����STRIP���Ȳ���SIDС��; �����߶�/��ȣ� ����LID��
function order = getITEMorder(Item,whichSortItemOrder)

%��SID����: SID������˳������,���С����ǰ��
szRow = cellfun(@(x)size(x,1), Item.SID);
if (max(szRow)~=min(szRow)),  error('ͬһITEM��Ӧ���ж��SID');  end %ͬһItemӦ��ֻ��һ��SID
SIDorder = cell2mat(Item.SID);   %ֱ��cell2matת��; %ITEM��SID 1-n��˳�򷵻� 

%��LID����: LID��ָ��˳��, ����SID����ȫ��һ��,�ٰ�LID��С��������,��ʵû������(��SID/LID����ͬһITEM),��󿴸߶� 
szRow = cellfun(@(x)size(x,1), Item.LID);
if (max(szRow)~=min(szRow)),  error('ͬһITEM��Ӧ���ж��SID');  end %ͬһItemӦ��ֻ��һ��LID
LIDorder = cell2mat(Item.LID);   %ֱ��cell2matת��; %ITEM��SID 1-n��˳�򷵻� 

% LIDorder = ones(1,length(SIDorder));
% *********** ����isNonMixed
tmpItem = [SIDorder; LIDorder; Item.LWH; Item.isNonMixed];  % tmpItem = [Item.SID; Item.LID; Item.LWH];  % tmpItem = [ Item.LWH];
        % [~,order] = sortrows(tmpItem',[1, 4, 2, 5],{'ascend','descend','descend','descend'}); 

% ����isNonMixed
    % [~,order] = sortrows(tmpItem',[1,6,  4, 3, 2, 5],{'ascend','descend','descend','descend','descend','descend'});  
    % �����ж�Item�ȳ���,��߶�����. % �����nbcol2����:
    % 1: SID ; 2: isNonMixed; 3: Longth/Height; 4: Height 5:Width; 6: LID; 
    % [~,order] = sortrows(tmpItem',[1,6, 4, 5, 3, 2],{'ascend','descend','descend','descend','descend','descend'});  
% Ŀǰ˳�� : 1: SID ; 2: isNonMixed; 3: Longth/Height; 4:Width; 5: LID; 6: Height
[~,order] = sortrows(tmpItem',[1,6, 4, 3, 2, 5 ],{'ascend','descend','descend','descend','descend','descend'});  
% [~,order] = sortrows(tmpItem',[1,6, 2, 4, 3, 5 ],{'ascend','descend','descend','descend','descend','descend'});  

% *********** ������isNonMixed
% % tmpItem = [SIDorder; LIDorder; Item.LWH; ];  % tmpItem = [Item.SID; Item.LID; Item.LWH];  % tmpItem = [ Item.LWH];
% % % [~,order] = sortrows(tmpItem',[1, 4, 2, 5],{'ascend','descend','descend','descend'}); 
% % [~,order] = sortrows(tmpItem',[1, 4, 3, 2, 5],{'ascend','descend','descend','descend','descend'});  

%��ITEM����/�����/���߶� ���������� ������ͬIDLU���ֿ�
% ����LUID 2: ȷ����ʹ������ȫ��ͬ ��LUID��ͬ�� Ҳ�����һ��
if ~isrow(order), order=order'; end
end

%% ����2: getThisLevel
function [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem,  Strip, p)
% ��ͬwhichStripH��,��ù�ͬ��thisLevel
if p.whichStripH == 1 % 1 bestfit 2 firstfit 3 nextfit
    % ����Rotation���ӱ���
    if sItem.isRota(iItem) == 1 %��Item������ת
        %             if ParaArray.whichRotation == 1
        % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem) |  ...
            Strip.LW(1, 1 : iLevel) >= sItem.LWH(2,iItem));
    else %��Item��������ת
        % ���������µ�ѡ��find����㹻�Ķ��level,����������Сʣ��ˮƽ��ȵ�
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem));
    end
    if isempty(flag)
        iLevel = iLevel + 1;% �����Ȳ����㣬��level����
        [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem,  Strip, p);
    else
        % ��ȡthisLevel: Ψһ��FF������⵽thisLevel�ļ��㣨ѡ��������������С�ģ�
        tmpLevels = Strip.LW(1,1:iLevel);   %��ȡ�����Ѱ��Ż��°��ŵ�level��ʣ��ˮƽƽ�������tmpLevels
        tepMinLeftWdith = min(tmpLevels(flag));                      %�ҳ�tepAvailableLevelArray�п����ɱ�iITem��ʣ��ˮƽ��ȵ���Сֵ������ΪtepMinLeftWdith
        thisLevel = find(tmpLevels==tepMinLeftWdith);            %�ҳ�����Сֵ��Ӧ���Ǹ�/Щlevel
        if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');          end
        if length(thisLevel)>1
            thisLevel = thisLevel(1);
        end
    end
elseif p.whichStripH == 2 % firstfit  % firstfit���ܲ���ֱ������bestfit�Ĵ���?
    % ����Rotation���ӱ���
    if sItem.isRota(iItem) == 1 %��Item������ת
        % �ҵ�����rotation�µ�level:��һ�ڷŷ���ɷ����iItem��level
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem) |  ...
            Strip.LW(1, 1 : iLevel) >= sItem.LWH(2,iItem));
    else
        % ���������µ�ѡ��find����㹻�Ķ��level,�������ڵ�һ�������� Ψһ������thisLevel�Ļ�ȡ
        flag = find(Strip.LW(1, 1 : iLevel) >= sItem.LWH(1,iItem));
    end
    if isempty(flag)
        iLevel = iLevel + 1;% �����Ȳ����㣬��level����
        [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem, Strip, p);
    else
        thisLevel = flag(1);
        if ~all(ismember(thisLevel,flag)),     error('Not all thisLevel belongs to flag ');       end
    end
elseif p.whichStripH == 3 % nextfit
    % ����Rotation���ӱ���
    if sItem.isRota(iItem) == 1 %��Item������ת % nextfit�²���ֱ������bestfit�Ĵ���
        % �ж���ǰlevel�Ƿ��������һ�ڷŷ���ɷ����iItem flaged: �����ݱ�ʾ���ԣ����򲻿���
        flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
        %                 flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) |  ...
        %                                       Strip.LW(1,iLevel) >= sItem.LWH(2,iItem));
    else
        % ��ͬ�����µ�ѡ�������ǰitem�Ŀ�<=��ǰstrip�ĵ�ǰlevel�Ŀ�
        flaged = find(Strip.LW(1,iLevel) >= sItem.LWH(1,iItem) );
    end
    if  isempty(flaged)  %ע����֮ǰ~flag������
        iLevel = iLevel + 1;% �����Ȳ����㣬��level����
        [thisLevel,iLevel] = getThisLevel( iItem, iLevel, sItem, Strip, p) ;
    else
        if  isempty(flaged) ,   error(' �����ܵĴ��� ');      end
        thisLevel = iLevel; % ��ǰlevelһ���ŵ���
    end
    %             end
    
end
end
