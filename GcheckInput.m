%% GCHECKINPUT ��Ҫ����:��������daת��+��������da�˶�
%% Form
%    [d] = GcheckInput(d,ParaArray)
%% Descripttion
%    �������ݼ���
%
%% Inputs
%    d                      (1,1)    �ṹ�壺�û���������
%    ParaArray          (1,1)    �ṹ�壺�㷨���Բ���
%
%% Outputs
%    d                     (1,1)
%

function [d] = GcheckInput(d)
    % Check if INPUT is legal with minor modification
     [d.LU, d.Veh] = initCheck(d.LU,d.Veh); %DO More time

    % Deep Check if INPUT is logical legal
    deepCheck(d.LU,d.Veh);
end

function deepCheck(LU,Veh)
    % �߶�Լ�� ����
    if min(Veh.LWH(3, :)) < max(LU.LWH(3, : )),     error('����: �������̸߶� ���� �����Ϳ��ø߶����� \n');     end

    % ����Լ�� ���� (���̳������)
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


    %% Ƕ�׺���
%     %% ����ת������
%     function initCheck()
%         % �ж�С��0
%         fLU = structfun(@(x) any( x(:) < 0), d.LU, 'UniformOutput', false);
%         fVeh = structfun(@(x) any( x(:) < 0), d.LU, 'UniformOutput', false);
%         if any(cell2mat(struct2cell(fLU))) || any(cell2mat(struct2cell(fVeh)))
%             error('�������ݴ���С��0����'); end
%         
%         nVeh = length(d.Veh.ID);
%         nLU = length(d.LU.ID);
%         % �ж�Veh�����Ƿ�������; �������Ƿ�һֱ
%         fields  = fieldnames(d.Veh);
%         for idx = 1:length(fields)
%             aField = d.Veh.(fields{idx});
%             if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���'); end
%             if isvector(aField)
%                 if ~isrow(aField)
%                     d.Veh.(fields{idx}) = aField';    end
%                 if length(aField) ~= nVeh
%                     error('Veh����������ڴ���');   end
%             else % ismatrix
%                 if size(aField,1) == nVeh && size(aField,2) ~= nVeh
%                     d.Veh.(fields{idx}) = aField';
%                 if size(aField,2)  ~= nVeh
%                     error('Veh����������ڴ���');   end                    
%                 end
%             end                
%         end
%         % �ж�LU�����Ƿ�������; �������Ƿ�һֱ
%         fields  = fieldnames(d.LU);
%         for idx = 1:length(fields)
%             aField = d.LU.(fields{idx});
%             if ~ismatrix(aField) || ~isnumeric(aField), error('�������ݴ��ڴ���'); end
%             if isvector(aField)
%                 if ~isrow(aField)
%                     d.LU.(fields{idx}) = aField';     end
%                 if length(aField) ~= nLU
%                     error('LU����������ڴ���');   end
%             else % ismatrix
%                 if size(aField,1) == nLU && size(aField,2) ~= nLU
%                     d.LU.(fields{idx}) = aField';    end
%                 if size(aField,2)  ~= nLU
%                     error('LU����������ڴ���');   end
%             end
%         end
%     end




    %% LUID����ת��: d.LU.ID ת��Ϊ��1��ʼ������� ������������ID��Ϣ
%     function idExchange()
% 
%         
% %         uniLUID = unique(d.LU.ID);
% %         nLUid = length(uniLUID);
% %         tmpLUID=d.LU.ID; %�м����
% %         for i=1:nLUid
% %             tmpLUID(d.LU.ID(:)==uniLUID(i)) = i;
% %         end
% %         d.LU.ID=tmpLUID;
%     end