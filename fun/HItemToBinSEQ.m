%% V2
function [LU,Item,Strip] = HItemToBinSEQ(LU,Item,Strip)

nItem = size(Item.LWH,2); %具体使用的Item的数量
nLU= size(LU.LWH,2); %具体使用的LU的数量

% 新增变量初始化
Item.Item_Bin=zeros(2,nItem);
LU.LU_Bin=zeros(2,nLU);

[LU.leftA,LU.bottomA,LU.rightA,LU.topA] = deal(zeros(1,nLU)); %LU的四个坐标

iiStrip = 0;
% 循环每个bin   
for iBin=1:max(Strip.Strip_Bin(1,:))
    tmpItemSeq=1;  %每个bin内顺序从1开始
    tmpLUSeq=1;
    tmpStrip = []; %计算CoordItemBin长(高)度使用 %       tmpLWStrip = Strip.LW(:,Strip.Strip_Bin(1,:)==iBin); %此bin内剩余宽度和长(高)度
    
    nbStrip=numel(find(Strip.Strip_Bin(1,:)==iBin));    %nbStrip：该iBin内strip数量   
    
    % 循环bin内每个strip
    for iStrip=1:nbStrip %从安放的第一个Strip开始逐步放置
        iiStrip = iiStrip+1;
        tmpLUSeqinStrip=1;
        [~,thisStrip] = find(Strip.Strip_Bin(1,:)==iBin & Strip.Strip_Bin(2,:)==iStrip ); %此bin内第iStrip个strip        

 
        if ~isscalar(thisStrip),     error('意外错误');     end
        
        nbItem=numel(find(Item.Item_Strip(1,:)==thisStrip));
        % 循环strip每个item,  从第一个顺序开始逐步获取ITEM在BIN内的坐标, 依据ITEM_STRIP的顺序
        for iItem=1:nbItem %同一strip内 循环取值iItem
            
            [~,thisItem] = find(Item.Item_Strip(1,:)==thisStrip & Item.Item_Strip(2,:)==iItem );    %此Strip内第iItem个item
            if ~isscalar(thisItem), error('意外错误');  end
            
            % 更新itemBeBinMatrix 555
            Item.Item_Bin(1,thisItem)=iBin;  % 本thisItem所在iBin的序号
            Item.Item_Bin(2,thisItem)=tmpItemSeq;  % 本thisItem在本iBin内的顺序 从1开始
            tmpItemSeq=tmpItemSeq+1;    
                    
            % 增加对LU的更新
            tmpLU = [];             %计算CoordLUBin高度使用 % tmpLWLU = LU.LWH(:,LU.LU_Item(1,:)==thisItem); 
            nbLU=numel(find(LU.LU_Item(1,:)==thisItem));     
            
            % 循环item每个LU, 从第一个顺序开始逐步获取LU在BIN内的坐标
            for iLU=1:nbLU                
                [~,thisLU] = find(LU.LU_Item(1,:)==thisItem & LU.LU_Item(2,:)==iLU);
                if ~isscalar(thisLU), error('意外错误');  end
                
                                            % 更新LU_Strip % NOTE: 在plot3DBPP时 出错, 发现此处赋值无用, so 删除该语句
                                            %                 LU.LU_Strip(1,thisLU)=iiStrip;
                                            %                 LU.LU_Strip(2,thisLU)=tmpLUSeqinStrip;

                tmpLUSeqinStrip=tmpLUSeqinStrip+1;
                                            % 更新LURotaed 555  % LURotaed(1,thisLU)=Item.Rotaed(1,thisItem);
                % 更新LU_Bin 555
                LU.LU_Bin(1,thisLU)=iBin;
                LU.LU_Bin(2,thisLU)=tmpLUSeq;
                tmpLUSeq=tmpLUSeq+1;
                
          tmpLU = [tmpLU thisLU];
                
           
            end
        end  % EOF ITEM
        tmpStrip = [tmpStrip thisStrip];
    end  %EOF STRIP
end  %EOF BIN




end % EOF function
