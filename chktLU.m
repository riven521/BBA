%% ���ܣ������Ƿ��������������أ��Ƿ�Z�߶���ItemSEQһ�£��Ƿ�ITEMID��XY����һ��;
function chktLU(t) % t�ض����������͵�table��struct
%% 0 �ṹ��ת��Ϊtable Ԥ����
if isstruct(t)
    t = struct2table(structfun(@(x) x', t,'UniformOutput',false));  end

% t.Properties.VariableNames
                                                                            % tLU1 = t(:,{'ID','LWH','LID','SID','PID','Weight',,'LU_Bin','LU_Item','isShuaiWei','Rotaed'});
%% 1 CHECK Weight �������� ͨ����ͬITEMID ��X Y �����ж� Ӧ�ǽ��׼
% 1.1 �����˲�LU_Item �粻���ڻ�ȡ��ITEMID 
if any(strcmp('LU_Item', t.Properties.VariableNames))
        t.ITEMID = t.LU_Item(:,1);
        t.ITEMSEQ = t.LU_Item(:,2); end

                                                                                                                            % if any(strcmp('ITEMID', t.Properties.VariableNames))
                                                                                                                            %     if any(t.LU_Item(:,1) ~= t.ITEMID) || any(t.LU_Item(:,2) ~= t.ITEMSEQ) %�Ƿ���ڼ����� �������ϰ汾, ������LU_ItemΪ��׼;
                                                                                                                            %             error('LU_Item��ITEMID��ITEMSEQ��ͬ'); end
if any(strcmp('CoordLUBin', t.Properties.VariableNames))
    t.X = t.CoordLUBin(:,1);    t.Y = t.CoordLUBin(:,2);   t.Z = t.CoordLUBin(:,3);   end
    
%% 2 ��ͬITEMID�µ�CHEK
if ~any(strcmp('LU_Item', t.Properties.VariableNames))
    return;
else
uniItemID = unique(t.ITEMID(:));
% 2.1 �������CoordLUBin����,������ϵ�ж���ITEMID�Ĳ���(1:XYֵ(�����Ƿ���ͬ); 2:������Zֵ)
for iItem = 1:length(uniItemID) %��ITEM����ѭ��
    flagIdx = t.ITEMID==uniItemID(iItem); %�Ե�һItemIdȥ�߼�ֵ;
    if any(strcmp('Y', t.Properties.VariableNames))
        vX = t{flagIdx,'X'};
        vY = t{flagIdx,'Y'};
        if any(vX ~= vX(1)) || any(vY ~= vY(1))
            error('��ͬITEM,��X��Y�����λ'); end
        
        v = t(flagIdx,{'Z','Weight','ITEMSEQ'});
    else
        v = t(flagIdx,{'Weight','ITEMSEQ'});
    end
    
    v = sortrows(v,'ITEMSEQ');
    
    if any(strcmp('Y', t.Properties.VariableNames))
        % ����������ǵݼ���Z�߶Ȳ��ǵ���
        if ~issorted(v.Z,'ascend') || ~issorted(v.Weight,'descend')
            issorted(v.Z,'ascend');  issorted(v.Weight,'descend');
            error('��ͬITEM,���������ǵݼ���Z�߶Ȳ��ǵ���'); end
    else
        if ~issorted(v.Weight,'descend')
            uniItemID(iItem)
            v;
            error('��ͬITEM,���������ǵݼ�'); end
        end
end
end



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