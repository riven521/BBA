%% ���ܣ������Ƿ��������������أ��Ƿ�Z�߶���ItemSEQһ�£��Ƿ�ITEMID��XY����һ��;
function chktLU(T) % t�ض����������͵�table��struct
%% 0 �ṹ��ת��Ϊtable Ԥ����
if isstruct(T), T = getTableLU(T); end

%% 1 CHECK LU (TODO ������ڰ�initCheck�ŵ��˴�,��chk����)

%% 2 CHECK ITEM (1:Weight �������� 
if ~any(strcmp('LU_Item', T.Properties.VariableNames)) %��û��,������δ�Ѷ�,��ֱ�ӷ���
    return;
end

% 2.1 CHECK �Ѷ������� ��ϵ�Ƿ���ȷ LU_Item
chkA_B(T.LU_Item')

uniItem = unique(T.ITEMID);

% 2.1 CHECK �Ƿ�����LW��ID��һ��
if any(strcmp('L', T.Properties.VariableNames)) && any(strcmp('ID', T.Properties.VariableNames))
    
        for iItem = 1:length(uniItem)
        flagLUIdx = T.ITEMID==uniItem(iItem);
        
        subT = T(flagLUIdx,{'L','W','ID'});        
        
        if ~isscalar(unique(subT.L)) || ~isscalar(unique(subT.W)) || ~isscalar(unique(subT.ID))
            error('��ͬITEM,��L��W��ID��һ��'); end        
        
        end
    
end

% 2.2 CHECK �Ѷ����� ������������
for iItem = 1:length(uniItem) %��ITEM����ѭ��
    flagLUIdx = T.ITEMID==uniItem(iItem); %�Ե�һItemIdȥ�߼�ֵ;
    
    subT = T(flagLUIdx,{'Weight','ITEMSEQ'});
    
    if  isUpDownWeight(subT) %����������
        error('��ͬITEM,���������ǵݼ�');  % T.isWeightUpDown(flagLUIdx) = 1;
    end
    
    %     subT = sortrows(subT,'ITEMSEQ');    
    %     if ~issorted(subT.Weight,'descend')            
    %             error('��ͬITEM,���������ǵݼ�');  % T.isWeightUpDown(flagLUIdx) = 1;
    %     end
end

% 2.3 CHECK �Ѷ����� �Ƿ����� X Y Z ��ȷ 
% 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYZֵ)
if any(strcmp('Y', T.Properties.VariableNames))
    
    for iItem = 1:length(uniItem)
        flagLUIdx = T.ITEMID==uniItem(iItem);
        
        subT = T(flagLUIdx,{'X','Y','Z','ITEMSEQ'});
        
        subT = sortrows(subT,'ITEMSEQ');
        
        if ~isscalar(unique(subT.X)) || ~isscalar(unique(subT.Y))
            error('��ͬITEM,��X��Y�����λ'); end
        
        if ~issorted(subT.Z,'ascend')
            error('��ͬITEM,���߶Ȳ��ǵ���');
        end
        
    end
end



%% 3 LWH1(��margin) LWH2������margin�� �� OLWH����Rotaed��������mragin�� ����
    % v1 : �ϰ汾,���ڿɱ������V2ȡ��
    OLWH = T.OLWH';    
    Rotaed = T.Rotaed';
    margin= T.margin';
    
    OLWHtmp = OLWH;
    OLWH(1,Rotaed) = OLWHtmp(2,Rotaed);
    OLWH(2,Rotaed) = OLWHtmp(1,Rotaed);
    
    LWH1 = T.LWH';
    LWH2 = LWHunbuffer(LWH1, margin);
    
    if ~isequal(OLWH,LWH1) && ~isequal(OLWH,LWH2)
        error('t��LWH�ڱ仯ǰ���д�����Ҫ����');
    end
    
    % v2: Ҳ�Ǻ���LWH�Ƿ���ȷ��
    LWHOLD = T.OLWH';
    LWHNEW = LWHunbuffer(T.LWH', T.margin', T.Rotaed');
    if ~isequal(LWHOLD,LWHNEW) 
        error('t��LWH�ڱ仯ǰ���д�����Ҫ����');
    end
end
  


%% �ֲ�����
function chkA_B(array)
if size(array,1) ~= 2, error('�������'); end

a = array(1,:);
b = array(2,:);
uniA = unique(a);

if ~issorted(uniA,'strictascend') || uniA(1)~= 1
    error('��һ�зǴ�1��ʼ�ϸ����');
end

for i=1:length(uniA)
    uniB = sort(b(a == uniA(i) ));
    if isempty(uniB) || ~issorted(uniB,'strictascend') || uniB(1)~= 1
        error('�ڶ��зǴ�1��ʼ�ϸ����');
    end
end

end



%% v2:����repair��check
% % function chktLU(t) % t�ض����������͵�table��struct
% % %% 0 �ṹ��ת��Ϊtable Ԥ����
% % if isstruct(t)
% %     t = struct2table(structfun(@(x) x', t,'UniformOutput',false));  end
% % 
% % % t.Properties.VariableNames
% %                                                                             % tLU1 = t(:,{'ID','LWH','LID','SID','PID','Weight',,'LU_Bin','LU_Item','isShuaiWei','Rotaed'});
% % %% 1 CHECK Weight �������� ͨ����ͬITEMID ��X Y �����ж� Ӧ�ǽ��׼
% % % 1.1 �����˲�LU_Item �粻���ڻ�ȡ��ITEMID 
% % if any(strcmp('LU_Item', t.Properties.VariableNames))
% %         t.ITEMID = t.LU_Item(:,1);
% %         t.ITEMSEQ = t.LU_Item(:,2); end
% % 
% %                                                                                                                             % if any(strcmp('ITEMID', t.Properties.VariableNames))
% %                                                                                                                             %     if any(t.LU_Item(:,1) ~= t.ITEMID) || any(t.LU_Item(:,2) ~= t.ITEMSEQ) %�Ƿ���ڼ����� �������ϰ汾, ������LU_ItemΪ��׼;
% %                                                                                                                             %             error('LU_Item��ITEMID��ITEMSEQ��ͬ'); end
% % if any(strcmp('CoordLUBin', t.Properties.VariableNames))
% %     t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
% %     
% % %% 2 ��ͬITEMID�µ�CHEK
% % if ~any(strcmp('LU_Item', t.Properties.VariableNames))
% %     return;
% % else
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYֵ(�����Ƿ���ͬ); 2:������Zֵ)
% % for iItem = 1:length(uniItemID) %��ITEM����ѭ��
% %     flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('��ͬITEM,��X��Y�����λ'); end
% %         
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
% %     else
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});
% %     end
% %     
% %     v = sortrows(v,'ITEMSEQ');
% %     
% %     if any(strcmp('Y', t.Properties.VariableNames))
% %         % ����������ǵݼ���Z�߶Ȳ��ǵ���
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend');  issorted(v.Weight,'descend');
% %             error('��ͬITEM,���������ǵݼ���Z�߶Ȳ��ǵ���'); end
% %     else
% %         if ~issorted(v.Weight,'descend')
% %             uniItemID(iItem)
% %             v;
% %             error('��ͬITEM,���������ǵݼ�'); end
% %         end
% % end
% % end
% % 
% % %% 3 LWH1(��margin) LWH2������margin�� �� OLWH����Rotaed��������mragin�� ����
% %     OLWH = t.OLWH';
% %     Rotaed = t.Rotaed';
% %     margin= t.margin';
% %     
% %     
% %     OLWHtmp = OLWH;
% %     OLWH(1,Rotaed) = OLWHtmp(2,Rotaed);
% %     OLWH(2,Rotaed) = OLWHtmp(1,Rotaed);
% %     
% %     LWH1 = t.LWH';
% %     LWH2 = LWHunbuffer(LWH1, margin);
% %     
% %     if ~isequal(OLWH,LWH1) && ~isequal(OLWH,LWH2)
% %         error('t��LWH�ڱ仯ǰ���д�����Ҫ����');
% %     end
% %     
%% V1 ���жϺ�ѭ��;
% % %% 2 ��ͬITEMID�µ�CHEK
% % uniItemID = unique(t.ITEMID(:));
% % % 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYֵ(�����Ƿ���ͬ); 2:������Zֵ)
% % if any(strcmp('CoordLUBin', t.Properties.VariableNames))
% %     if ~any(strcmp('Y', t.Properties.VariableNames))
% %         t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
% % if any(strcmp('Y', t.Properties.VariableNames))
% %     for iItem = 1:length(uniItemID) %��ITEM����ѭ��
% %         % a �������ж�
% %         flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
% %         vX = t{flagIdx,'X'};
% %         vY = t{flagIdx,'Y'};
% %         if any(vX ~= vX(1)) || any(vY ~= vY(1))
% %             error('��ͬITEM,��X��Y�����λ'); end
% %         % b �������������ж�
% %         v = t(flagIdx,{'Z','Weight','ITEMSEQ'});   v = sortrows(v,'ITEMSEQ');
% %         % ����������ǵݼ���Z�߶Ȳ��ǵ���
% %         if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
% %             issorted(v.Z,'ascend')
% %             issorted(v.Weight,'descend')
% %             error('��ͬITEM,���������ǵݼ���Z�߶Ȳ��ǵ���'); end
% %     end
% % end
% % % 2.2 �������CoordLUBin������,������ж�ITEMID�Ĳ���(2:����)
% % else 
% %     for iItem = 1:length(uniItemID)
% %         flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
% %         v = t(flagIdx,{'Weight','ITEMSEQ'});   v = sortrows(v,'ITEMSEQ');
% %         % ����������ǵݼ���Z�߶Ȳ��ǵ���
% %         if ~issorted(v.Weight,'descend')
% %             issorted(v.Weight,'descend')
% %             error('��ͬITEM,���������ǵݼ�'); end
% %     end
% % end