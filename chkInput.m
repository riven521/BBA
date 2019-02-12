function [d] = chkInput(d)
% GCHECKINPUT ��Ҫ����:��������daת��(�����������ݣ�+��������da�˶�
% ����1��initCheck
%   1��ȷ��LU��Veh��Ϊ��ֵ�;���
%   2���Եߵ���Veh����ת��
%   3���Բ�ͬSID������ͬID�Ľ���ID�޸ģ�REMOVE)
%   4���Եߵ���LU����ת��
% ����2��deepCheck
%   1��LU�߶�Լ��
%   2��LU����Լ��
%   3��LU�ɶѶ��һ����Լ��
     
     %  Check if INPUT is legal with minor modification
     [d.LU, d.Veh] = initCheck(d.LU,d.Veh); 

    % Deep Check if INPUT is logical legal
    deepCheck(d.LU,d.Veh);
end

% �ֲ�����
function [LU,Veh] = initCheck(LU,Veh)
% initCheck if INPUT is legal AND ת�õߵ��� INPUT 
% 1��ȷ��LU��Veh��Ϊ��ֵ�;���
% 2���Եߵ���Veh����ת��
% 3���Բ�ͬSID������ͬID�Ľ���ID�޸�
% 4���Եߵ���LU����ת��

nVeh = length(Veh.Weight);
nLU = length(LU.Weight);

%% *************** 1 Lu��Veh���ۺ��жϣ�ȷ��Ϊ��ֵ�;��� ***************
% 1: ������Ϊ�Ǹ�/��ά���飬��������������������Ϳ������double������
try
    structfun(@(x) validateattributes(x,{'double'},{'nonnegative','2d'}), LU, 'UniformOutput', false);
    structfun(@(x) validateattributes(x,{'double'},{'nonnegative','2d'}), Veh, 'UniformOutput', false);
catch
    fLU = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
    fVeh = structfun(@(x) any( x(:) < 0), LU, 'UniformOutput', false);
    if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
        error('�������ݴ���С��0����'); end
end

%% *************** 2 Veh�����ж� ***************
% �ж�Veh��ÿ��field�����ж�������/���� �Ƿ�һ�£��Ƿ�������ΪnVeh
fields  = fieldnames(Veh);
for idx = 1:length(fields)
    aField = Veh.(fields{idx});
    try
        % V2 �汾�ж�
        % �����ж�
        if ~strcmp('LWH',fields{idx}) % && ~strcmp('buff',fields{idx}) Ŀǰû��buff����
            if ~isrow(aField),  Veh.(fields{idx}) = aField';     end %����������ݷ�������,�޸�
            validateattributes(Veh.(fields{idx}),{'double'},{'nonnegative','vector','ncols', nVeh});
        % �����ж�
        else
            if size(aField,1) == nVeh && size(aField,2) ~= nVeh %������������==������ �� ����~=������
                Veh.(fields{idx}) = aField'; %����ת��
            end
            
            % ��������ж�������������Ϊ3�����������ж���
            if size(aField,1) == nVeh && size(aField,2) == nVeh
                
                if nVeh==3
                    Veh.(fields{idx}) = aField';
                    warning('���ڴ���3�����������,δ�����ж�,������ж����������Ƿ���ȷ');
                    fprintf('�����ĵĳ����Ϊ: %d \n ',Veh.(fields{idx}));
                    tmpVeh = Veh.(fields{idx});
                    if ~(tmpVeh(1)*1.5 <= tmpVeh(2)) %�����һ�������Ŀ�ȵ�1.5�����ڳ����ĳ���,����Ϊ����߷ŷ���
                        Veh.(fields{idx}) = aField;
                    end
                    
                else
                    warning('����3���������Ƿ���,�������ڷǳ����ʱ����');
                end
                
            end
            
            validateattributes(Veh.(fields{idx}),{'double'},{'nonnegative','2d','ncols', nVeh});
            
        end
    catch    
        error('try error in Veh�����ж�');
        % V1 �汾�ж�
        % �ж�1: ��������ֵ��
        if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���(����ֵ�;���)'); end
        % �ж�2: �����粻��nVeh,����ת��(�����������,���ڿ����Ǿ������͵�LWH��buff�ֿ��ж�)
        if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('buff',fields{idx}) %��������������Ҳ��Ǻ�������ֵ(�˴����жϳ��������LWH)
            if ~isrow(aField) %����������ݷ�������,�޸�
                Veh.(fields{idx}) = aField';         end
        else % ismatrix %�˴��жϳ���������Ƿ�д�� LWH buff��
            if size(aField,1) == nVeh && size(aField,2) ~= nVeh %������������==������ �� ����~=������
                Veh.(fields{idx}) = aField'; %����ת��
            end
            if size(aField,1) == nVeh && size(aField,2) == nVeh

                if nVeh==3
                    Veh.(fields{idx}) = aField';
                    warning('���ڴ���3�����������,δ�����ж�,������ж����������Ƿ���ȷ');
                    fprintf('�����ĵĳ����Ϊ: %d \n ',Veh.(fields{idx}));
                    tmpVeh = Veh.(fields{idx});
                    if ~(tmpVeh(1)*1.5 <= tmpVeh(2)) %�����һ�������Ŀ�ȵ�1.5�����ڳ����ĳ���,����Ϊ����߷ŷ���
                        Veh.(fields{idx}) = aField;
                    end             

                else
                    warning('����3���������Ƿ���,�������ڷǳ����ʱ����');  
                end

            end        
        end
        % �ж�3: ת�������������nVeh
        if size(Veh.(fields{idx}),2)  ~= nVeh %5555 ��Ҫ�ж�
            error('Veh���ݴ������������ڳ��������,�������2');
        end
    
    end
end

%% *************** 3 LU��ID�Ų����ڲ�ͬSID���ظ�(Milkrun����) - �ı�ID�� �ŵ�Ԥ������Gpreproc *************** 
     % 2 �����ͬID��PID/EID/LID�ȣ����£���ӦSID��Ҫ���벻ͬ������ͬ�ı�ID�ţ�ֱ����������ͬ��ID�ڲ�ͬSID��;
%      if isrepeated(LU.ID,LU.SID)
%         warning('��������ID���ظ�, ��Ҫ����'); 
%         LU.ID = reviseID(LU.ID,LU.SID);
%      end

    % V1: ������ɸ��º�ID������ID����ͬ
%     done = false;
%     while ~done
%     flag = 0;
%     uniID = unique(LU.ID);                  %��ͬ��ID��
%     for iID = 1:length(uniID)
%         fID = LU.ID==uniID(iID);
%         uniSID = unique(LU.SID(fID));
%         if length(uniSID) > 1                %������ڣ��ı�ID�ţ�ID+SID�ţ�
%             for iSID = 1:length(uniSID)
%                 fSID = LU.SID==uniSID(iSID);
%                 f = fSID & fID;
%                 LU.ID(f) = LU.ID(f) + LU.SID(f);               
%             end
%             flag = 1;                                   % ���flag==1���������޸ģ�����Ҫwhileѭ����ֱ����ѭ��
%         end
%     end
%     if flag==0, done=true; end
%     end
    
%% *************** 4 LU�����ж� ***************
% �ж�LU�����Ƿ�������; �������Ƿ�һ��; ��һ�½�����ת
fields  = fieldnames(LU);
for idx = 1:length(fields)
    aField = LU.(fields{idx});
      try
        % V2 �汾�ж�
        % �����ж�
        if ~strcmp('LWH',fields{idx}) && ~strcmp('maxL',fields{idx}) && ~strcmp('margin',fields{idx}) %Ŀǰû��buff����
            if ~isrow(aField),  LU.(fields{idx}) = aField';    end %����������ݷ�������,�޸�
            validateattributes(LU.(fields{idx}),{'double','logical'},{'nonnegative','vector','ncols', nLU});
        % �����ж�
        else
            if size(aField,1) == nLU && size(aField,2) ~= nLU
                LU.(fields{idx}) = aField';   end
            
            % ��������ж�������������Ϊ3�����������ж���
            if size(aField,1) == nLU && size(aField,2) == nLU
                
                if nLU==3
                    LU.(fields{idx}) = aField';
                    warning('���ڴ���3�����̵����,���������');
                    fprintf('���̵ĵĳ����Ϊ: %d \n ',LU.(fields{idx}));
                    uniID = unique(LU.ID);
                    % ����ID - ��Ӧ����LW/isRota/yID/maxL  etc... ���������ת��
                    for iLU = 1:length(uniID)
                        tmpM = LU.(fields{idx});
                        tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %��ȡ��ͳ�,��Ҫ��
                        if isempty(tmpVeh),   error('��Ӧ�ó��ֵĴ���');        end
                        if numel(unique(tmpVeh,'rows')) > 2
                            LU.(fields{idx}) = aField;
                            warning('ԭ���̳���ߵȽ�����ת�� \n');
                            break;
                        end
                    end
                else
                    warning('����3���������Ƿ���,�������ڷǳ����ʱ����');
                end
                
            end
            
            validateattributes(LU.(fields{idx}),{'double'},{'nonnegative','2d','ncols', nLU});
            
        end
    catch
         error('try error in LU�����ж�');
        % V1 �汾�ж�
        % �ж�1: ��������ֵ��
        if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���(����ֵ�;���)'); end
        % �ж�2: �����粻��nLU,����ת��    (�����������,���ڿ����Ǿ������͵�LWH��buff��margin�ֿ��ж�)
        if isvector(aField) && ~strcmp('LWH',fields{idx}) && ~strcmp('maxL',fields{idx}) && ~strcmp('margin',fields{idx})%&& ~strcmp('buff',fields{idx} ) %��������������Ҳ��Ǻ�������ֵ
            if ~isrow(aField)
                LU.(fields{idx}) = aField';     end
        else % ismatrix %�˴��ж����̳�����Ƿ�д��
            if size(aField,1) == nLU && size(aField,2) ~= nLU
                LU.(fields{idx}) = aField';   end
            
            if size(aField,1) == nLU && size(aField,2) == nLU
                
                if nLU==3
                    LU.(fields{idx}) = aField';
                    warning('���ڴ���3�����̵����,���������');
                    fprintf('���̵ĵĳ����Ϊ: %d \n ',LU.(fields{idx}));
                    uniID = unique(LU.ID);
                    % ����ID - ��Ӧ����LW/isRota/yID/maxL  etc... ���������ת��
                    for iLU = 1:length(uniID)
                        tmpM = LU.(fields{idx});
                        tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %��ȡ��ͳ�,��Ҫ��
                        if isempty(tmpVeh),   error('��Ӧ�ó��ֵĴ���');        end
                        if numel(unique(tmpVeh,'rows')) > 2
                            LU.(fields{idx}) = aField;
                            warning('ԭ���̳���ߵȽ�����ת�� \n');
                            break;
                        end
                    end
                else
                    warning('����3���������Ƿ���,�������ڷǳ����ʱ����');
                end
                
            end
        end
        % �ж�3: ת�������������nLU
        if size(LU.(fields{idx}),2)  ~= nLU
            error('LU���ݴ����������������������,���ڲ���ԭ�µĴ���2');
        end
    end
end 

end

function deepCheck(LU,Veh)
    % LU�߶�Լ�� ����
    if min(Veh.LWH(3, :)) < max(LU.LWH(3, : )),     error('����: �������̸߶� ���� �����Ϳ��ø߶����� \n');     end

    % LU����Լ�� ���� (���̳������)
    flag1 = LU.isRota == 1;  %����rota��LU���
    flag0 = LU.isRota == 0;  %������rota��LU���
    
    % ���Բ�������ת��LU
    if min(Veh.LWH(1, :)) < max(LU.LWH(1,flag0)), error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
    if min(Veh.LWH(2, :)) < max(LU.LWH(2,flag0)), error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
    
    % ����������ת��LU
    if ~all(~flag1) && min(Veh.LWH(1:2)) < max(min(LU.LWH(1:2, flag1)))  %flag1����ȫ0 �������������
        error('����: ����������̱ߵ����ֵ ���� �����ͳ���� \n');
    end
    if ~all(~flag1) && max(Veh.LWH(1:2)) < max(max(LU.LWH(1:2, flag1 ))) %flag1����ȫ0 �������������
        error('����: ����������ߵ����ֵ ���� �����ͳ���� \n');
    end

    %% V2 : ����LU����������ID������ͬ ���Ҳ�ͬ��Ӧ�̣� ���䳤��ͬ����
% %     uniSID = unique(LU.SID);
% %     uniID = unique(LU.ID);
% %     % ����ID - ��Ӧ����/isRota/yID/maxL  ����555
% %     for iSID = 1:length(uniSID)
% %         fSID = LU.SID==uniSID(iSID);
% %         for iLU = 1:length(uniID)
% %             fID = LU.ID==uniID(iLU);
% %             f = fSID & fID;
% %             tmpLULW = LU.LWH(1:2, f)';                       if isempty(tmpLULW),  error('��Ӧ�ó��ֵĴ���');  end %��ȡ��ͳ�,��Ҫ��           
% %             if numel(unique(tmpLULW,'rows')) > 2,  error('����: ��������ID������ͬ ���䳤��ͬ���� \n');    end
% %             
% %             tmpLULW = LU.isRota(:,f)';                           if isempty(tmpLULW),      error('��Ӧ�ó��ֵĴ���');       end
% %             if numel(unique(tmpLULW)) > 2,      error('����: ��������ID������ͬ ����������ת���Ͳ�ͬ���� \n');    end
% %             
% %             %         tmp = LU.maxL(:,LU.ID==uniID(iLU));
% %             %         if isempty(tmp),      error('��Ӧ�ó��ֵĴ���');       end
% %             %         if numel(unique(tmp)) > 2,  error('����: ��������ID������ͬ ����������maxL���Ͳ�ͬ���� \n');    end
% %             %         tmp = LU.yID(:,LU.ID==uniID(iLU));
% %             %         if isempty(tmp),      error('��Ӧ�ó��ֵĴ���');       end
% %             %         if numel(unique(tmp)) > 2,  error('����: ��������ID������ͬ ����yID���Ͳ�ͬ���� \n');    end
% %         end
% %     end
    

    %% V1 : ����LU����������ID������ͬ ���䳤��ͬ���� (ͨ���޸�ID,ȷ����ͬSID�µ�ID�ǲ�ͬ��)
    % (ȷ������ͬSID��Ӧ���µ����̲���ƴ��)
    uniID = unique(LU.ID);
    % ����ID - ��Ӧ����/isRota/yID/maxL  ����555
    for iLU = 1:length(uniID)
        tmp = LU.LWH(1:2,LU.ID==uniID(iLU))';   %��ȡ��ͳ�,��Ҫ�� 
        if isempty(tmp),      error('��Ӧ�ó��ֵĴ���');       end
                        %x=[LU.LWH;LU.ID]
                %         x=[LU.LWH;LU.ID;LU.SID]
                %         f = LU.ID==uniID(iLU)
                %         y=x(:,f)'
                %         unique(y,'rows')
        if numel(unique(tmp,'rows')) > 2,  error('����: ��������ID������ͬ ���䳤��ͬ���� \n');    end        
        tmp = LU.isRota(:,LU.ID==uniID(iLU));
        if isempty(tmp),      error('��Ӧ�ó��ֵĴ���');       end
        if numel(unique(tmp)) > 2,  error('����: ��������ID������ͬ ����������ת���Ͳ�ͬ���� \n');    end
%         tmp = LU.maxL(:,LU.ID==uniID(iLU));
%         if isempty(tmp),      error('��Ӧ�ó��ֵĴ���');       end
%         if numel(unique(tmp)) > 2,  error('����: ��������ID������ͬ ����������maxL���Ͳ�ͬ���� \n');    end
%         tmp = LU.yID(:,LU.ID==uniID(iLU));
%         if isempty(tmp),      error('��Ӧ�ó��ֵĴ���');       end
%         if numel(unique(tmp)) > 2,  error('����: ��������ID������ͬ ����yID���Ͳ�ͬ���� \n');    end
    end
    

%% INPUT CHECK: ��������ͬLID����ͬSID����ͬEID��ȴ����ƴ�أ���ID��ͬ���ƺ�û���������
% %     tlu = array2table([LU.SID;LU.EID;LU.LID]')
% %     tluid = array2table([LU.ID]')
% %     [a,b]=findgroups(tlu);
% %     x = splitapply(@(x) length(unique(x)),tluid.Var1,a)
% % %     Num = splitapply(@(w) length(unique(w)),T2.DUEDATA,GID) %ȡֵ��������
% %    


end



%                         tmpVeh =tmpM(1:2,LU.ID==uniID(iLU))';    %��ȡ��ͳ�,��Ҫ��

    
%%
                                    %         fLU = structfun(@isrow, rmfield(d.LU,'LWH'), 'UniformOutput', false);
                                    %         fVeh = structfun(@isrow, rmfield(d.Veh,{'LWH','buff'}), 'UniformOutput', false);
                                    %         idx = find(cell2mat( struct2cell(fVeh) )== 0)

                                    % �����ж�����: d.Veh
                                        %         if numel(d.Veh.LWH(:))~=numel(d.Veh.BUFF(:)) || numel(d.Veh.LWH(:)) ~= 3 ||...
                                        %                 ~isscalar(d.Veh.Weight)
                                        %             error('������da.BinArray����,����Ԥ�� \n');
                                        %         end
                                            % �����ж�����: d.LU
                                    %         if numel(d.LU.ID(:))~=numel(d.LU.Weight(:)) || numel(d.LU.BUFF(:)) ~= 2 || ...
                                    %                 numel(d.LU.ID(:))*3~=numel(d.LU.LWH(:)) || ~ismatrix(d.LU.LWH)    
                                    %             error('������da.LUArray����,����Ԥ�� \n');

% if any(dBin(1) < LU.LWH(1,flag0)),  error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
% if any(dBin(2) < LU.LWH(2,flag0)),  error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
% ����Լ�� ���� (���̳������) TOBE DELE
%         if ParaArray.whichRotation == 1
%             if min(dBin(1:2)) < max(min(dLU(1:2,:))) || max(dBin(1:2)) < max(max(dLU(1:2,:)))
%                 error('����: �������̳���� ���� �����ͳ���� \n');
%             end
%         else %��׼rotation
%             if any(dBin(1) < dLU(1,:)),  error('����: �������̿�� ���� �����Ϳ��ÿ������ \n');    end
%             if any(dBin(2) < dLU(2,:)),  error('����: �������̳��� ���� �����Ϳ��ó������� \n');    end
%         end