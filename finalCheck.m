%% ���Table����
function Tseq = finalCheck(T,D)
% CoordLUBin  LWH  LU_VehType    BINID    BINSEQ    SID    LID    ITEMID    PID    ShowSEQ    Weight    tblorder
% ID    LID            LWH             isRota           margin           PID    SID    Weight
%% 0 oD ��ԭʼ���������Ľ��
oD = D(T.tblorder,:);
% % %% 1.1 LWH CHECK
% % T.LWH(T.Rotaed,1) = oD.LWH(T.Rotaed,1);
% % T.LWH(T.Rotaed,2) = oD.LWH(T.Rotaed,2);
% % if any(T.LWH ~= oD.LWH), any(T.LWH ~= oD.LWH); 
% %     T.LWH
% %     oD.LWH
% %     error('����1');  end
% % %% 1.2 Weight CHECK
% % if any(T.Weight ~= oD.Weight), any(T.Weight ~= oD.Weight);  
% %     T.Weight 
% %     oD.Weight 
% %     error('����2');  end
%% 2 ����Tseq���� LID PID SID 
T.SID = oD.SID;
T.PID = oD.PID;
T.LID = oD.LID;
Tseq = T(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});

T1=TLUsorted(:,'CoordLUBin');
T2=TLUsorted(:,'LWH');
T3= TLUsorted(:,{'LU_VehType','BINID','BINSEQ','ITEMID','ShowSEQ','Rotaed','tblorder','Weight'});

%%
% T1 = T(:,{'LWH'});
% T2 = D(:,{'LWH'});
% ismember(T.LWH(:,:),D.LWH(T.tblorder,:)) %������
% ismember(T1,T2)

% if any(T.LWH ~= D.LWH),    error('����');  end
% T123.OPID == T123.OPID
% D.OPID == T123.OPID
% T1 = T123(:,{'LID','OPID','OSID','Weight'});
% T2 = D(:,{'LID','OPID','OSID','Weight'});
% ismember(T.LWH(:,:),TIN.LWH(T.tblorder,:))
% T1.PID
% T2.PID
% if ~all(ismember(T2,T1))
%     ismember(T2,T1)
%     setdiff(T2,T1)
%     setdiff(T1,T2)
% %     setdiff(A(:,vars),B(:,vars))
%     error('����'); 
% end