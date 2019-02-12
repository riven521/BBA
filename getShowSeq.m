function [TLU] = getShowSeq(d)

if ~istable(d)
    TLU = getTableLU(d);
else
    TLU=d;
end

%% 1 Table���� : TLUsorted %����ǰ����������Ѿ����� BINSEQһ���̶��ģ�ֻ��ShowSEQ�ֲ������
[TLUsorted,tblorder]= sortrows(TLU,{'BINID','BINSEQ','SID','LID','PID','ITEMID','ITEMSEQ','H','LU_VehType'},...
                                         {'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend','ascend'});
%% 2 ��������չʾ˳��ShowSEQ����ֵ��TLUsorted
        % ����LUShowSeq : ���һ��: REAL����չʾ˳��(��˦β��)
        % Ŀǰ����LU_Bin(1,:)�ֳ�, SID�ֹ�Ӧ��, LID���������� ��������  TODO ��������������Ҫ�жϲ�������� % [2 4 5]        
        tmpThreeRows = TLUsorted(:,{'BINID','SID','LID'});
        LUShowSeq=zeros(height(tmpThreeRows),1);
        LUShowSeq(1)=1;
        for i =2:length(LUShowSeq)
            if tmpThreeRows{i,'BINID'}~=tmpThreeRows{i-1,'BINID'} %����ֳ�,�����ʼ��1
                LUShowSeq(i) = 1;
            else %�������SID��LID����ͬ,���������
                if isequal(tmpThreeRows(i,{'SID','LID'}),tmpThreeRows(i-1,{'SID','LID'}))
                    LUShowSeq(i) = LUShowSeq(i-1) ;
                else
                    LUShowSeq(i) = LUShowSeq(i-1)+1;
                end
            end            
        end        

%% ����TLU
       %         TLUsorted.ShowSEQ = LUShowSeq;        %         TLUsorted.tblorder = tblorder;
       [~,x] = sort(tblorder);
       TLU.ShowSEQ = LUShowSeq(x);  %        TLU.ShowSEQ = LUShowSeq(tblorder);
       TLU.tblorder = tblorder;
       
%% 4 ����Ϊ���ص����ݱ�
% T1=TLUsorted(:,'CoordLUBin');
% T2=TLUsorted(:,'LWH');
% T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','ITEMID','ShowSEQ','Rotaed','tblorder','Weight'}); %ɾ���̶�ֵ,��TLUIN��ȡ
% % T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});

end




%% V2 : V3 ׼��ɾ������Ҫ��ע��
% % function [TLU] = getShowSeq(d)
% % 
% % if ~istable(d)
% %     TLU = getTableLU(d);
% % else
% %     TLU=d;
% % end
% % 
% % %% 1 Table���� : TLUsorted %����ǰ����������Ѿ����� BINSEQһ���̶��ģ�ֻ��ShowSEQ�ֲ������
% % [TLUsorted,tblorder]= sortrows(TLU,{'BINID','BINSEQ','SID','LID','PID','ITEMID','ITEMSEQ','H','LU_VehType'},...
% %                                          {'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend','ascend'});
% % %% 2 ��������չʾ˳��ShowSEQ����ֵ��TLUsorted
% %         % ����LUShowSeq : ���һ��: REAL����չʾ˳��(��˦β��)
% %         % Ŀǰ����LU_Bin(1,:)�ֳ�, SID�ֹ�Ӧ��, LID���������� ��������  TODO ��������������Ҫ�жϲ�������� % [2 4 5]        
% %         tmpThreeRows = TLUsorted(:,{'BINID','SID','LID'});
% %         LUShowSeq=zeros(height(tmpThreeRows),1);
% %         LUShowSeq(1)=1;
% %         for i =2:length(LUShowSeq)
% %             if tmpThreeRows{i,'BINID'}~=tmpThreeRows{i-1,'BINID'} %����ֳ�,�����ʼ��1
% %                 LUShowSeq(i) = 1;
% %             else %�������SID��LID����ͬ,���������
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
% % %% 3 ����Ϊ���ص����ݱ�
% % % T1=TLUsorted(:,'CoordLUBin');
% % % T2=TLUsorted(:,'LWH');
% % % T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','ITEMID','ShowSEQ','Rotaed','tblorder','Weight'}); %ɾ���̶�ֵ,��TLUIN��ȡ
% % % % T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});
% % 
% % %% 4 �ƺ����� ���ڿ�ɾ
% % % % s1= table2array(T1)';
% % % % s2= table2array(T2)';
% % % % s3= table2array(T3,'ToScalar',true)';
% % % % % if ~isequal(output_CoordLUBin,s1) || ~isequal(output_LU_LWH,s2) || ~isequal(output_LU_Seq(1:7,:),s3(1:7,:))
% % % % %     error('V1 and V2����ͬ');
% % % % % end
% % % % output_CoordLUBin=s1;
% % % % output_LU_LWH=s2;
% % % % output_LU_Seq=s3;
% % 
% % end















%% V1 % ���ز���1��2��3 ���ڽṹ��
% % function [output_CoordLUBin,output_LU_LWH,output_LU_Seq] = getReturnBBA(daMax)
% % %% 1 ����������(ԭʼ˳��) ���3������
% % 
% % % ����1 - LU��Bin�ڵ�����
% % % V2:  LU margin��ʽ
% % output_CoordLUBin = daMax.LU.CoordLUBin;
% % 
% %     % V1:  LU buff ��϶��ʽ
% %     % daMax.LU.CoordLUBinWithBuff = daMax.LU.CoordLUBin + daMax.LU.buff./2;
% %     % output_CoordLUBin=daMax.LU.CoordLUBinWithBuff; %output_CoordLUBin��DOUBLE����: Lu��xyzֵ TTTTTTTTTT
% % 
% % % ����2 - LU�ĳ����(��ת��)
% % % LWH�Ѿ�Ϊ��С�����Ӧmargin���ʵ�����ݱ���
% % % ������V3 - LU margin��ʽ
% % output_LU_LWH = daMax.LU.LWH; %output_LU_LWH��DOUBLE LU�ĳ���ߣ���ת��ʵ��ֵ��
% % 
% %         % ������V2
% %         %  ���Ӽ�϶-�޶�LWHΪ��С�����ӦBuffer���ʵ�����ݱ���
% %         %  daMax.LU.LWHOriRota = daMax.LU.LWH - daMax.LU.buff;
% %         %  output_LU_LWH=daMax.LU.LWHOriRota;  %output_LU_LWH��DOUBLE LU�ĳ���ߣ���ת��ʵ��ֵ��
% %         % ������V1
% %         %         daMax.LU.LWHRota = daMax.LU.LWHRota - daMax.LU.BUFF;
% %         %         Res3_LWHRota=daMax.LU.LWHRota;  %Res3_LWHRota��DOUBLE LU�ĳ���ߣ���ת��
% % 
% % %% ����3 - GET sortedDaMax: ��С���ȵ�ԪLUչʾ�ľۺ�
% % % 3.1 �������������� (��δ����LU_VehType)
% %     LU_Item=daMax.LU.LU_Item;
% %     LID=daMax.LU.LID;                   %LU�Ѷ���LUID, ������˳����LID % LID=daMax.LU.ID;
% %     PID=daMax.LU.PID;
% %     SID=daMax.LU.SID;
% %     hLU=daMax.LU.LWH(3,:);
% %     LU_Bin = daMax.LU.LU_Bin;   %Ψһ���е�
% %     LU_VehType=daMax.LU.LU_VehType;
% % 
% %     if ~isfield(daMax.LU, 'isShuaiWei')
% %         LShuaiWei = zeros(length(LU_VehType),1)';
% %     else
% %         LShuaiWei = daMax.LU.isShuaiWei;
% %     end
% %     % LPingPu = daMax.LU.isPingPu; �ݲ�����
% %      
% %     
% % sortedDaMax = [LU_Item; LID; PID; SID; hLU; LU_Bin;LU_VehType]; % 2 1 1 1 1 2 1
% % % �����Ҫ��LUID���㲿���󰴶Ѷ�չʾ, ȡͬһBIN��, ͬһSID, ͬһLUID ->> ͬһ PID, ͬһLU_ITEM
% % % 1 BIN 2 BINSEQ 3 SID 4 LID -> 5 PID 6 ITEM 7 ITEMSEQ 8 LUHEIGHT 7==8 9 LU_VehType
% % sortedDaMax = [LU_Bin(1,:); LU_Bin(2,:); SID; LID; PID; ...
% %                            LU_Item(1,:); LU_Item(2,:); hLU; LU_VehType; LShuaiWei];
% % 
% %         % ��������˳�� tmpSeq:
% %         % �����Ҫ��LUID�ȶѶ�չʾ,���㲿��չʾ, ȡͬһBIN��, ͬһSID, ͬһLUID ->> ͬһ LU_ITEM��ͬһPID
% %         % 1 BIN 2 BINSEQ 3 SID 4 LID -> 5 ITEM 6 ITEMSEQ 7 PID 8 LUHEIGHT 6==8 
% %         %         tmpSeq =[7,8,5,3,1,2,4,6];
% %         % �����Ҫ��LUID���㲿���󰴶Ѷ�չʾ, ȡͬһBIN��, ͬһSID, ͬһLUID ->> ͬһ PID, ͬһLU_ITEM
% %         % 1 BIN 2 BINSEQ 3 SID 4 LID -> 5 PID 6 ITEM 7 ITEMSEQ 8 LUHEIGHT 7==8 9 LU_VehType
% %             % tmpSeq =[7,8,5,3,4,1,2,6];
% %             % [~,order] = sortrows(sortedDaMax',tmpSeq,{'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend'});
% % 
% % % V2 ������
% % tmpSeq = [1:9];
% % [~,order] = sortrows(sortedDaMax',tmpSeq,{'ascend','ascend','ascend','ascend','ascend','ascend','ascend','descend','ascend'});
% % 
% % % FINAL return's results; ��order˳�򷵻�
% % output_LU_LWH = output_LU_LWH(:,order);
% % output_CoordLUBin = output_CoordLUBin(:,order);
% % 
% % sortedDaMax =sortedDaMax(:,order);
% % 
% % %% ����3 - GET output_LU_Seq: 
% % % 3.2 ��������չʾ˳��
% %         % ����LUShowSeq : ���һ��: REAL����չʾ˳��(��˦β��)
% %         % Ŀǰ����LU_Bin(1,:)�ֳ�, SID�ֹ�Ӧ��, LID���������� ��������  TODO ��������������Ҫ�жϲ�������� % [2 4 5]
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
% % % 3.2 ���������н��չʾ��Щ�м���˳��:
% % % 9 LU_VehType 1 BIN 2 BINSEQ 3 SID 4 LID 6 ITEM(LU_Item(1,:)) 5 PID
% % % ��1���������ڳ��ͺ�    ��2���������ڳ���� ��3�����̳��ڰ���˳�� ��4������SID��Ӧ�̱��
% % % ��5������ID�ͺ�LID ��6�����̶Ѷ����ITEM ��7�������㲿�����PID ������8: չʾ˳��
% % 
% % % V3 ��Ҫ���صľ���ֵ
% % order
% %     LU_VehType = sortedDaMax(9,:);   %��1���������ڳ��ͺ�          LU_VehType;
% %     LUBINID = sortedDaMax(1,:);             %��2���������ڳ����      LU_Bin(1,:)
% %     LUBINSeq = sortedDaMax(2,:);        %��3�����̳��ڰ���˳��     LU_Bin(2,:)
% %     LUSID = sortedDaMax(3,:);        %��4������SID��Ӧ�̱��           SID
% %     LULID = sortedDaMax(4,:);        %  ��5������ID�ͺ�LID                LID
% %     LUItemID = sortedDaMax(6,:);        % ��6�����̶Ѷ����ITEM      LU_Item(1,:)
% %     LUPID = sortedDaMax(5,:);        % ��7�������㲿�����PID         PID
% %     LUShowSeq                                 % ��8: ����չʾ˳��                  SHOWSEQ
% %     output_LU_Seq = [LU_VehType; LUBINID; LUBINSeq; LUSID; ...
% %                                 LULID; LUItemID; LUPID; LUShowSeq];
% %     
% %                 %  V1 1 LU_VehType 2 LU_Bin(1) 3 LU_Bin(2) 4 SID 5 LID 6 LU_Item(1) 7 PID
% %                 % tmpShow =[9,7,8,5,3,1,4];  %����9:�����������ͺ� ����3���к�     % tmpShow =[7,8,5,3,1,2,4,6];
% %                 % V2  % ���7��: ��1: LU_VehType ��2 LU_Bin(1) ��3 LU_Bin(2) ��4 SID ��5 LID ��6
% %                 %  LU_Item(1) ��7 PID ��8 ���˳��(����Ҫ)
% %                 % 1-9: [LU_Bin(1,:); LU_Bin(2,:); SID; LID; PID; LU_Item(1,:);
% %                 % LU_Item(2,:); hLU; LU_VehType]; ����10 LShuaiWei
% %                 % tmpShow =[9,1,2,3,4,6,5,10];
% %                 %output_LU_Seq =sortedDaMax(tmpShow,:);
% %                 %output_LU_Seq = [sortedDaMax(1,:);LUShowSeq];
% % 
% % end
