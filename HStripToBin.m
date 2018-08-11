function [Strip,Bin]= HStripToBin(Strip,Veh,LU,p)
% 重要函数:Strip放入Bin中 %  行数:长宽高(row);  列数:托盘数量(coloum);
% Input ---  Strip/Veh:  
% Output --- Strip: 
% Output --- Bin: 

%% 初始化
sz = size(Strip.LW);
nStrip = sz(2);

wVeh  = Veh.LWH(1,1); 
lVeh  = Veh.LWH(2,1); 

%% Strip排序 555 (如何确保相同托盘类型的相邻摆放？？？ TODO FIX ME)
    % 获取Strip的顺序(重点是Strip高度递减排序（但经常遇到strip高度一样的）) %
    [Strip.striporder] = getStriporder(Strip);  % Strip两排序方式 高度/长度递减
    % 获取按order排序后的Strip: sStrip    
    if isSameCol(Strip)
        sStrip = structfun(@(x) x(:,Strip.striporder),Strip,'UniformOutput',false);
    else
        error('不能使用structfun');
    end
    
   
%% LU->Item->Strip->Bin转换 
% 获取stripBeBinMatrixSort: 每个排序后strip在哪个bin内  以及顺序
% 获取LWBin:  新生成的Bin的剩余长宽
Bin.LW = zeros(2,nStrip);    %初始化bin: dim1-bin宽度剩余 ; dim2-bin长(高)度(555剩余）;
Bin.LW(1,:) = wVeh;
Bin.LW(2,:) = lVeh;
Bin.Weight = zeros(1,nStrip); % 初始赋值

    % 初始化多行nItem列
    Bin.LID = zeros(numel(unique(LU.ID)),nStrip);
    Bin.PID = zeros(numel(unique(LU.PID)),nStrip);
    Bin.SID = zeros(numel(unique(LU.SID)),nStrip);
    Bin.UID = zeros(numel(unique(LU.UID)),nStrip);
    
tmpBin_Strip = zeros(1,nStrip);    % 每个Bin内的Strip数量 后期不用
% sStrip新增
sStrip.Strip_Bin = zeros(2,nStrip); % dim1:序号 strip在某个bin dim2:进入顺序 555

% LWStripSort = sStrip.LW; %从sorted获得
% StripWeightSort = sStrip.Weight;

% 55 获取thisBin - 当前strip要放入的bin序号
% 循环往bin中安置strip,即固定strip,变化选择不同bin(thisBin)
% 注释：获取 FLAG        可放下当前iStrip的至少一个bin的集合 
% 注释：获取 thisBin   从FLAG中找到按规则的那个thisBin, 并执行 insert函数

iStrip=1; iBin=1;
while 1
    if iStrip > nStrip, break; end
    
    [thisBin,iBin] = getThisBin(iBin, iStrip, sStrip, Veh, Bin, p);    % 获取Bin号
    
    [Bin,sStrip.Strip_Bin,tmpBin_Strip] = insertStripToBin(iStrip, thisBin, sStrip, Bin, sStrip.Strip_Bin, tmpBin_Strip);
        
    iStrip = iStrip + 1;
end

% plot2DBin();

% 后处理 并赋值到da
% 获取Strip_Bin: 每个strip在哪个bin内  以及顺序
                % Strip.Strip_Bin( : , Strip.striporder) = sStrip.Strip_Bin;
% Strip内部更新,sStrip依据order变化回来
if isSameCol(sStrip)
    Strip = reorderStruct(Strip.striporder, sStrip);
else
    error('不能使用structfun');
end

% 获取Bin: 去除未使用的Bin 注意Bin结构体的变化
if isSameCol(Bin)
    Bin = structfun(@(x) x( : , Bin.Weight(1,:)>0 ), Bin, 'UniformOutput', false);
else
    error('不能使用structfun');
end


   
%% 测试script
% 输出主要结果:获得每个item包含的 原始 LU序号
printscript();
    

%% 嵌套函数


    function printscript()
        % 输出主要结果:获得从1开始每个bin包含的数据
        % Strip.stripBeBinMatrix
        for iBin = 1:max(Strip.Strip_Bin(1,:))
            [~,idx] = find(Strip.Strip_Bin(1,:)==iBin); %本iBin下的strip索引号
            idxSeq = Strip.Strip_Bin(2,idx); %本iBin内strip放入顺序Seq
            fprintf('bin 的宽+长为: ' );
            fprintf(' %d  ', Veh.LWH(:,1));
            fprintf('\n');
            fprintf('bin %d 的剩余宽+剩余长为:  ',iBin);
            fprintf('( %d ) ',Bin.LW(:,iBin));
            fprintf('\n');
            fprintf('bin %d 包含 original strip 索引号{顺序}(长宽)为  \n  ',iBin);
            fprintf('%d ',idx);fprintf('\n');
            fprintf('{%d} ',idxSeq);fprintf('\n');
            fprintf('( %d ) ', Strip.LW(1:2,idx));fprintf('\n');
            fprintf('\n');
        end
    end

    % 未完成函数 TODO
    function plot2DBin()
    % 初始化
            % 初始化
        Bin.LW
        sStrip.LW
        tmpBin_Strip
        sStrip.Strip_Bin
        sStrip
            %% 初始化
        nThisItem = size(d.Item.LWH,2);
        nIDType = unique(d.Item.LID);
        nColors = hsv(length(nIDType)); %不同类型LU赋予不同颜色        
%         tmpUniqueBin = unique(Veh.LWH(1:2,:)','rows')';
        %         wBin = tmpUniqueBin(1);
%         hBin = tmpUniqueBin(2);     
        wBin = Veh.LWH(1,1);
        hBin = Veh.LWH(2,1);
   
    
        nUsedBin = sum(sStrip.Strip_Bin(2,:)>0);

%         %% 画图
        % 1 画个画布 宽度为nUsedBin+1个bin宽 长（高）度为bin高
        figure();
        DrawRectangle([wBin*(nUsedBin+1)/2 hBin/2 wBin*(nUsedBin+1) hBin 0],'--');
        hold on;
        % 2 逐个bin 画图
        iterWidth=0;    %每个bin在前1个bin的右侧 此为增加变量
        for iBin = 1:nUsedBin
            % 找出当前iBin的物品索引
            idxDrawStrip = find(sStrip.Strip_Bin(1,:)==iBin);
            % 。。。 由于没有Strip在bin内的Coord，此函数暂停
        end
        % 逐个strip画图
        
    end

end


    function order = getStriporder(Strip)
%         tmpLWStrip = Strip.LW(1:2,:);
%         [~,order] = sort(tmpLWStrip(2,:),'descend');  %对strip进行排序,只需要它的顺序ord;按第nDim=2行排序（长/高度)

Strip.SID
zs = getOrderofID(Strip.SID)
Strip.LID
zl = getOrderofID(Strip.LID)
%  t = Strip.SID;
%  ss = sum(t);  %每个STRIP内包含的SID个数
%  
%   for i=1:size(t,1)
%      if sum(ss( find(t(i,:))  ) > 1) > 1
%          error('有同一个SID被2个及以上STRIP包括');
%      end
%  end
% if  any(sum(t)>2)
%     error('有同一个strip包括3个及以上各SID');
% end
% 
%  z = zeros(1,size(t,2)); 
%  k=1;
%  for i=1:numel(ss)
%     if ss(i) ==1
%         if i>1 && ss(i-1) ==1 && find(t(:, i)==1) ~= find(t(:, i-1)==1) %判断当前与前一个strip是否属于同样SID
%             k=k+1;
%         end
%         z(i) = k;
%     elseif ss(i) >1 %只要遇到STRIP包含2个及以上的STRIP时，更新顺序        
%         k=k+1;
%         z(i)=k;
%         k=k+1;
%     end
%  end

                     %  ss = sum(tSID);
                    %  torder = zeros(1,length(sorder));
                    %  k=1;
                    %  for i=1:length(sorder)
                    %      tt = tSID(sorder(i),:)
                    %      to
                    % %      tSID(sorder(i)) = [];
                    %      other = sorder;
                    %      other(sorder(i)) = [];
                    %      tti = tSID(sorder(other),:);
                    %      torder(tt==1) = k;
                    %      k = k+1;
                    %      f = tti==1 & tt==1 %本次SID有，但其它里面也有,排序为k+1
                    %      if any(f)
                    %          torder(f) = k; 
                    %          k = k+1;
                    %      end
                    %      1
                    %  end
                    % 
                    % Strip.SID
                    % [a,b,~]=find(Strip.SID==1)
                    % [a,b,~]=find(Strip.SID(:,:)==1)
                    % Strip.LID

        tmpSort = [zs; zl; Strip.LW(1:2,:); Strip.loadingrateLimit;Strip.loadingrate];
%         [~,order] = sortrows(tmpSort',[1],{'ascend'}); %对strip进行排序;按第nDim=2行排序（长/高度)，再看strip内部loadingrateLimit
        [~,order] = sortrows(tmpSort',[1,2],{'ascend','ascend'}); %对strip进行排序;按第nDim=2行排序（长/高度)，再看strip内部loadingrateLimit
% descend ascend
%         tmpSort = [zorder; Strip.LW(1:2,:); Strip.loadingrateLimit;Strip.loadingrate];
%         [~,order] = sortrows(tmpSort',[1 3 4 5 ],{'ascend','descend','descend','descend'}); %对strip进行排序;按第nDim=2行排序（长/高度)，再看strip内部loadingrateLimit
       
%         tmpLWH = [tmpIDItem; tmpLWHItem]; %额外增加ITEM的ID到第一行形成临时变量
%         [~,order] = sortrows(tmpLWH',[3 1 ],{'descend','descend'}); %按高度,ID(相同高度时)递减排序
%         tmpLWH = [tmpLWHItem; min(tmpLWHItem(1:2,:))]; %额外增加最短边到第三行形成临时变量tmpLWH
%         [~,order] = sortrows(tmpLWH',[3 2 1],{'descend','descend','descend'}); %way1 按最短边,高度,宽度递减排序
      
        if ~isrow(order), order=order'; end
    end

        function [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p)        
        if p.whichBinH == 1 % 1 bestfit
            % 条件: 寻找 bin的剩余高度 >= 本strip的高度 &且 bin的剩余重量 >= 本strip的重量 (集合中的最小值)
            flag = find(Bin.LW(2,1:iBin) >= sStrip.LW(2,iStrip)  & ...
                Veh.Weight(1) - Bin.Weight(1 : iBin) >= sStrip.Weight(iStrip) ); %
            if isempty(flag)
                iBin = iBin + 1; % 如果高度不满足，则bin升级
                [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p);                
            else
                tepBins = Bin.LW(2,1 : iBin); %获取所有已安排或新安排的bin的剩余高度向量tepBins
                tepMin = min(tepBins(flag)); % 555 check 找出bin中能放istrip且高度最小值tepMin（TODO 是否考虑重量？）
                thisBin = find(Bin.LW(2,1:iBin)==tepMin); %找到该值tepMin对应的那个/些bin序号
                if ~all(ismember(thisBin,flag)),      error('Not all thisBin belongs to flag ');        end
                if length(thisBin)>1
                    thisBin = thisBin(1);
                end
            end
        elseif p.whichBinH == 2 % 1 firstfit
            flag = find(Bin.LW(2,1:iBin) >= sStrip.LW(2,iStrip)  & ...
                Veh.Weight(1) - Bin.Weight(1 : iBin) >= sStrip.Weight(iStrip) );            
            if isempty(flag)
                iBin = iBin + 1; 
                [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p);    
            else
                thisBin = flag(1);
                if ~all(ismember(thisBin,flag)),     error('Not all thisBin belongs to flag ');       end
            end
        elseif p.whichBinH == 3 % 1 nextfit
            flaged = find(Bin.LW(2, iBin) >= sStrip.LW(2,iStrip) & ...
                Veh.Weight(1) - Bin.Weight(iBin) >= sStrip.Weight(iStrip) );            
            if  isempty(flaged)  %注意与之前~flag的区别
                iBin = iBin + 1; 
                [thisBin,iBin]  = getThisBin( iBin, iStrip, sStrip, Veh, Bin, p);    
            else
                if  isempty(flaged) ,   error(' 不可能的错误 ');      end
                thisBin = iBin; % 当前bin一定放的下
            end
        else
            error('错误参数设置');
        end
        end

        
    
    function [Bin,Strip_Bin,Bin_Strip] = insertStripToBin(iStrip, thisBin,sStrip,Bin,Strip_Bin,Bin_Strip)
%         binBeStripArray=binBeStripArray;stripBeBinMatrixSort=stripBeBinMatrixSort;Bin.LW=Bin.LW;

        % 1 更新Bin相关非Sort数据        
        %  1.1 更新strip归属bin的信息 ：stripBeBinMatrixSort
        Bin_Strip(thisBin) = Bin_Strip(thisBin) + 1; %本bin下第几次安置strip
        
        %  1.2 更新本bin的剩余长和剩余高：Bin.LW
        Bin.LW(1,thisBin) = min(Bin.LW(1,thisBin),sStrip.LW(1,iStrip)); %更新bin剩余宽度的最小值
        Bin.LW(2,thisBin) = Bin.LW(2,thisBin) - sStrip.LW(2,iStrip);    %更新bin剩余高度
            
        %  1.3 更新本bin对应的BinWeight: 
        Bin.Weight(thisBin) =  Bin.Weight(thisBin) + sStrip.Weight(iStrip);
        
        % 2 更新Strip相关Sort数据
        %  2.1 更新stripBeBinMatrixSort
        Strip_Bin(1,iStrip) = thisBin;
        Strip_Bin(2,iStrip) = Bin_Strip(thisBin);
   
            % 更新bIN中包含ID类与否
            Bin.LID(:,thisBin) = Bin.LID(:,thisBin) + sStrip.LID(:,iStrip); % 数值为出现次数
            Bin.LID(Bin.LID>0) = 1; % 数值改为出现与否
            Bin.PID(:,thisBin) = Bin.PID(:,thisBin) + sStrip.PID(:,iStrip); % 数值为出现次数
            Bin.PID(Bin.PID>0) = 1; % 数值改为出现与否
            Bin.SID(:,thisBin) = Bin.SID(:,thisBin) + sStrip.SID(:,iStrip); % 数值为出现次数
            Bin.SID(Bin.SID>0) = 1; % 数值改为出现与否
            Bin.UID(:,thisBin) = Bin.UID(:,thisBin) + sStrip.UID(:,iStrip); % 数值为出现次数
            Bin.UID(Bin.UID>0) = 1; % 数值改为出现与否
            
            
       %% 其余放到ItemToBin内计算
    end
