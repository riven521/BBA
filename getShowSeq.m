function [TLU] = getShowSeq(d)

if ~istable(d)
    TLU = getTableLU(d);
else
    TLU=d;
end

%% 1 Table排序 : TLUsorted %基本前两行排序就已经定了 BINSEQ一定固定的，只是ShowSEQ分步骤而已
[TLUsorted,tblorder]= sortrows(TLU,{'BINID','BINSEQ','SID','LID','PID','ITEMID','ITEMSEQ','H','LU_VehType'},...
                                         {'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend','ascend'});
%% 2 计算托盘展示顺序ShowSEQ并赋值到TLUsorted
        % 计算LUShowSeq : 最后一行: REAL托盘展示顺序(含甩尾等)
        % 目前按照LU_Bin(1,:)分车, SID分供应商, LID分托盘种类 三个区分  TODO 后期增加其它需要判断步骤的依据 % [2 4 5]        
        tmpThreeRows = TLUsorted(:,{'BINID','SID','LID'});
        LUShowSeq=zeros(height(tmpThreeRows),1);
        LUShowSeq(1)=1;
        for i =2:length(LUShowSeq)
            if tmpThreeRows{i,'BINID'}~=tmpThreeRows{i-1,'BINID'} %如果分车,步骤初始化1
                LUShowSeq(i) = 1;
            else %如果车内SID与LID均不同,步骤号增加
                if isequal(tmpThreeRows(i,{'SID','LID'}),tmpThreeRows(i-1,{'SID','LID'}))
                    LUShowSeq(i) = LUShowSeq(i-1) ;
                else
                    LUShowSeq(i) = LUShowSeq(i-1)+1;
                end
            end            
        end        

%% 返回TLU
       %         TLUsorted.ShowSEQ = LUShowSeq;        %         TLUsorted.tblorder = tblorder;
       [~,x] = sort(tblorder);
       TLU.ShowSEQ = LUShowSeq(x);  %        TLU.ShowSEQ = LUShowSeq(tblorder);
       TLU.tblorder = tblorder;
       
%% 4 精简为返回的数据表
% T1=TLUsorted(:,'CoordLUBin');
% T2=TLUsorted(:,'LWH');
% T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','ITEMID','ShowSEQ','Rotaed','tblorder','Weight'}); %删除固定值,从TLUIN获取
% % T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});

end




%% V2 : V3 准备删除不需要的注释
% % function [TLU] = getShowSeq(d)
% % 
% % if ~istable(d)
% %     TLU = getTableLU(d);
% % else
% %     TLU=d;
% % end
% % 
% % %% 1 Table排序 : TLUsorted %基本前两行排序就已经定了 BINSEQ一定固定的，只是ShowSEQ分步骤而已
% % [TLUsorted,tblorder]= sortrows(TLU,{'BINID','BINSEQ','SID','LID','PID','ITEMID','ITEMSEQ','H','LU_VehType'},...
% %                                          {'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend','ascend'});
% % %% 2 计算托盘展示顺序ShowSEQ并赋值到TLUsorted
% %         % 计算LUShowSeq : 最后一行: REAL托盘展示顺序(含甩尾等)
% %         % 目前按照LU_Bin(1,:)分车, SID分供应商, LID分托盘种类 三个区分  TODO 后期增加其它需要判断步骤的依据 % [2 4 5]        
% %         tmpThreeRows = TLUsorted(:,{'BINID','SID','LID'});
% %         LUShowSeq=zeros(height(tmpThreeRows),1);
% %         LUShowSeq(1)=1;
% %         for i =2:length(LUShowSeq)
% %             if tmpThreeRows{i,'BINID'}~=tmpThreeRows{i-1,'BINID'} %如果分车,步骤初始化1
% %                 LUShowSeq(i) = 1;
% %             else %如果车内SID与LID均不同,步骤号增加
% %                 if isequal(tmpThreeRows(i,{'SID','LID'}),tmpThreeRows(i-1,{'SID','LID'}))
% %                     LUShowSeq(i) = LUShowSeq(i-1) ;
% %                 else
% %                     LUShowSeq(i) = LUShowSeq(i-1)+1;
% %                 end
% %             end            
% %         end
% %         
% % % %          tmpThreeRows = TLUsorted(:,{'BINID','SID','LID'});
% % % %         LUShowSeq=zeros(height(tmpThreeRows),1);
% % % %         LUShowSeq(1)=1;
% % % %         for i =2:length(LUShowSeq)
% % % %             tmpThreeRows{i,'LID'}
% % % %             tmpThreeRows{i-1,'LID'}
% % % %             if tmpThreeRows{i,'BINID'}==tmpThreeRows{i-1,'BINID'} && ...
% % % %                     tmpThreeRows{i,'SID'}==tmpThreeRows{i-1,'SID'}  && ...
% % % %                     tmpThreeRows{i,'LID'}==tmpThreeRows{i-1,'LID'} 
% % % %                     LUShowSeq(i) = LUShowSeq(i-1) ;
% % % %                 else
% % % %                     LUShowSeq(i) = LUShowSeq(i-1)+1;
% % % %                 end
% % % %         end
% % 
% % %         TLUsorted.ShowSEQ = LUShowSeq;
% % %         TLUsorted.tblorder = tblorder;
% %        [~,x] = sort(tblorder);
% % %        TLU.ShowSEQ = LUShowSeq(tblorder);
% %        TLU.ShowSEQ = LUShowSeq(x);
% %        TLU.tblorder = tblorder;
% %        
% % %        TLUsorted.Weight
% % %        TLU.Weight
% % %        
% % % TLUsorted.ShowSEQ
% % %        TLU.ShowSEQ
% % %        1
% % 
% % %% 3 精简为返回的数据表
% % % T1=TLUsorted(:,'CoordLUBin');
% % % T2=TLUsorted(:,'LWH');
% % % T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','ITEMID','ShowSEQ','Rotaed','tblorder','Weight'}); %删除固定值,从TLUIN获取
% % % % T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});
% % 
% % %% 4 似乎无用 后期可删
% % % % s1= table2array(T1)';
% % % % s2= table2array(T2)';
% % % % s3= table2array(T3,'ToScalar',true)';
% % % % % if ~isequal(output_CoordLUBin,s1) || ~isequal(output_LU_LWH,s2) || ~isequal(output_LU_Seq(1:7,:),s3(1:7,:))
% % % % %     error('V1 and V2不相同');
% % % % % end
% % % % output_CoordLUBin=s1;
% % % % output_LU_LWH=s2;
% % % % output_LU_Seq=s3;
% % 
% % end















%% V1 % 返回参数1，2，3 基于结构体
% % function [output_CoordLUBin,output_LU_LWH,output_LU_Seq] = getReturnBBA(daMax)
% % %% 1 返回输出结果(原始顺序) 输出3个参数
% % 
% % % 参数1 - LU在Bin内的坐标
% % % V2:  LU margin方式
% % output_CoordLUBin = daMax.LU.CoordLUBin;
% % 
% %     % V1:  LU buff 间隙方式
% %     % daMax.LU.CoordLUBinWithBuff = daMax.LU.CoordLUBin + daMax.LU.buff./2;
% %     % output_CoordLUBin=daMax.LU.CoordLUBinWithBuff; %output_CoordLUBin：DOUBLE类型: Lu的xyz值 TTTTTTTTTT
% % 
% % % 参数2 - LU的长宽高(旋转后)
% % % LWH已经为减小长宽对应margin后的实际数据变量
% % % 以下是V3 - LU margin方式
% % output_LU_LWH = daMax.LU.LWH; %output_LU_LWH：DOUBLE LU的长宽高（旋转后：实际值）
% % 
% %         % 以下是V2
% %         %  增加间隙-修订LWH为减小长宽对应Buffer后的实际数据变量
% %         %  daMax.LU.LWHOriRota = daMax.LU.LWH - daMax.LU.buff;
% %         %  output_LU_LWH=daMax.LU.LWHOriRota;  %output_LU_LWH：DOUBLE LU的长宽高（旋转后：实际值）
% %         % 以下是V1
% %         %         daMax.LU.LWHRota = daMax.LU.LWHRota - daMax.LU.BUFF;
% %         %         Res3_LWHRota=daMax.LU.LWHRota;  %Res3_LWHRota：DOUBLE LU的长宽高（旋转后）
% % 
% % %% 参数3 - GET sortedDaMax: 最小粒度单元LU展示的聚合
% % % 3.1 参数三的列排序 (暂未考虑LU_VehType)
% %     LU_Item=daMax.LU.LU_Item;
% %     LID=daMax.LU.LID;                   %LU堆垛用LUID, 但返回顺序用LID % LID=daMax.LU.ID;
% %     PID=daMax.LU.PID;
% %     SID=daMax.LU.SID;
% %     hLU=daMax.LU.LWH(3,:);
% %     LU_Bin = daMax.LU.LU_Bin;   %唯一两行的
% %     LU_VehType=daMax.LU.LU_VehType;
% % 
% %     if ~isfield(daMax.LU, 'isShuaiWei')
% %         LShuaiWei = zeros(length(LU_VehType),1)';
% %     else
% %         LShuaiWei = daMax.LU.isShuaiWei;
% %     end
% %     % LPingPu = daMax.LU.isPingPu; 暂不考虑
% %      
% %     
% % sortedDaMax = [LU_Item; LID; PID; SID; hLU; LU_Bin;LU_VehType]; % 2 1 1 1 1 2 1
% % % 如果需要按LUID先零部件后按堆垛展示, 取同一BIN内, 同一SID, 同一LUID ->> 同一 PID, 同一LU_ITEM
% % % 1 BIN 2 BINSEQ 3 SID 4 LID -> 5 PID 6 ITEM 7 ITEMSEQ 8 LUHEIGHT 7==8 9 LU_VehType
% % sortedDaMax = [LU_Bin(1,:); LU_Bin(2,:); SID; LID; PID; ...
% %                            LU_Item(1,:); LU_Item(2,:); hLU; LU_VehType; LShuaiWei];
% % 
% %         % 排序优先顺序 tmpSeq:
% %         % 如果需要按LUID先堆垛展示,后零部件展示, 取同一BIN内, 同一SID, 同一LUID ->> 同一 LU_ITEM，同一PID
% %         % 1 BIN 2 BINSEQ 3 SID 4 LID -> 5 ITEM 6 ITEMSEQ 7 PID 8 LUHEIGHT 6==8 
% %         %         tmpSeq =[7,8,5,3,1,2,4,6];
% %         % 如果需要按LUID先零部件后按堆垛展示, 取同一BIN内, 同一SID, 同一LUID ->> 同一 PID, 同一LU_ITEM
% %         % 1 BIN 2 BINSEQ 3 SID 4 LID -> 5 PID 6 ITEM 7 ITEMSEQ 8 LUHEIGHT 7==8 9 LU_VehType
% %             % tmpSeq =[7,8,5,3,4,1,2,6];
% %             % [~,order] = sortrows(sortedDaMax',tmpSeq,{'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend'});
% % 
% % % V2 列排序
% % tmpSeq = [1:9];
% % [~,order] = sortrows(sortedDaMax',tmpSeq,{'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend','ascend'});
% % 
% % % FINAL return's results; 按order顺序返回
% % output_LU_LWH = output_LU_LWH(:,order);
% % output_CoordLUBin = output_CoordLUBin(:,order);
% % 
% % sortedDaMax =sortedDaMax(:,order);
% % 
% % %% 参数3 - GET output_LU_Seq: 
% % % 3.2 计算托盘展示顺序
% %         % 计算LUShowSeq : 最后一行: REAL托盘展示顺序(含甩尾等)
% %         % 目前按照LU_Bin(1,:)分车, SID分供应商, LID分托盘种类 三个区分  TODO 后期增加其它需要判断步骤的依据 % [2 4 5]
% %         tmpThreeRows = sortedDaMax([1,3,4],:);
% %         LUShowSeq=zeros(1,size(tmpThreeRows,2));
% %         LUShowSeq(1)=1;
% %         if length(LUShowSeq)>1
% %             for i =2:length(LUShowSeq)
% %                 if tmpThreeRows(1,i)==tmpThreeRows(1,i-1) && tmpThreeRows(2,i)==tmpThreeRows(2,i-1) && tmpThreeRows(3,i)==tmpThreeRows(3,i-1)
% %                     LUShowSeq(i) = LUShowSeq(i-1) ;
% %                 else
% %                     LUShowSeq(i) = LUShowSeq(i-1)+1;
% %                 end
% %             end
% %         end
% %         %         ThreeRows=[ThreeRows;LUShowSeq]
% % 
% % % 3.2 参数三的行结果展示哪些行及其顺序:
% % % 9 LU_VehType 1 BIN 2 BINSEQ 3 SID 4 LID 6 ITEM(LU_Item(1,:)) 5 PID
% % % 行1：托盘所在车型号    行2：托盘所在车序号 行3：托盘车内安置顺序 行4：托盘SID供应商编号
% % % 行5：托盘ID型号LID 行6：托盘堆垛序号ITEM 行7：托盘零部件编号PID 增加行8: 展示顺序
% % 
% % % V3 需要返回的具体值
% % order
% %     LU_VehType = sortedDaMax(9,:);   %行1：托盘所在车型号          LU_VehType;
% %     LUBINID = sortedDaMax(1,:);             %行2：托盘所在车序号      LU_Bin(1,:)
% %     LUBINSeq = sortedDaMax(2,:);        %行3：托盘车内安置顺序     LU_Bin(2,:)
% %     LUSID = sortedDaMax(3,:);        %行4：托盘SID供应商编号           SID
% %     LULID = sortedDaMax(4,:);        %  行5：托盘ID型号LID                LID
% %     LUItemID = sortedDaMax(6,:);        % 行6：托盘堆垛序号ITEM      LU_Item(1,:)
% %     LUPID = sortedDaMax(5,:);        % 行7：托盘零部件编号PID         PID
% %     LUShowSeq                                 % 行8: 托盘展示顺序                  SHOWSEQ
% %     output_LU_Seq = [LU_VehType; LUBINID; LUBINSeq; LUSID; ...
% %                                 LULID; LUItemID; LUPID; LUShowSeq];
% %     
% %                 %  V1 1 LU_VehType 2 LU_Bin(1) 3 LU_Bin(2) 4 SID 5 LID 6 LU_Item(1) 7 PID
% %                 % tmpShow =[9,7,8,5,3,1,4];  %增加9:托盘所出车型号 参数3的行号     % tmpShow =[7,8,5,3,1,2,4,6];
% %                 % V2  % 输出7行: 行1: LU_VehType 行2 LU_Bin(1) 行3 LU_Bin(2) 行4 SID 行5 LID 行6
% %                 %  LU_Item(1) 行7 PID 行8 输出顺序(最重要)
% %                 % 1-9: [LU_Bin(1,:); LU_Bin(2,:); SID; LID; PID; LU_Item(1,:);
% %                 % LU_Item(2,:); hLU; LU_VehType]; 增加10 LShuaiWei
% %                 % tmpShow =[9,1,2,3,4,6,5,10];
% %                 %output_LU_Seq =sortedDaMax(tmpShow,:);
% %                 %output_LU_Seq = [sortedDaMax(1,:);LUShowSeq];
% % 
% % end
