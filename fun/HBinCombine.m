function T = HBinCombine(do,flaggetSmallVeh,do1,flagTiledArray,do2Array,do3Array)
% HBinCombine ==> ��ԭʼ �ı䳵�� ƽ�� ��ϼ�϶ ��Ľ�����кϲ�

global ISlastVehType ISpingpu parGap

%% 1 ******************��ȡչʾ˳�� do���� -> T
    T = getTableLU(do); chktLU(T);

    % ��ͼ��ԭʼ��ƽ�̵�ͼ
    %     plotSolutionT(T,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));

%%
% ��1���������ڳ��ͺ�(���뻻)    ��2���������ڳ����(���,���ܻ�,���ʹ�) ��3�����̳��ڰ���˳��(���뻻) ��4������SID��Ӧ�̱��(�����,���ñ�??)
% ��5������ID�ͺ�LID(�����,���ñ�?) ��6�����̶Ѷ����ITEM(���,���ܻ�,���ʹ�,��;?) ��7�������㲿�����PID(�����,���ñ�?) ������8: չʾ˳��(���뻻)

%% 2 ****************** ��Գ��ͱ仯 do1���� ��ȡ�޶��� output ******************
    if ISlastVehType==1 && flaggetSmallVeh == 1 %���е������ҳ����滻�ɹ�
        T1 = getTableLU(do1);     chktLU(T1);

        % ��ȡ ���һ�����߼�ֵ 
        flaglastLUIdx = T.BINID == max(T.BINID);           % V1 V2�ȼ�: % lastVehIdx = max(T.BINID);  %lastVehIdx = max(T{:,'BINID'});     % flaglastLUIdx = T.BINID == lastVehIdx; % lastVehIdx = ibin

        %% ��Щ��仯?? fixeme ����ʹ��ʱ����ϸ����
        % ĳ��bin�ڵ���,��binIDһ������仯;
        % ��LID/Weight/LWHӦ�ò���仯; SID/PID��仯; ��ΪOPID OSID OID��ԭ��
        % ĳ��bin�ڵ���,��BINSEQ,CoordLUBin,LU_VehTypeһ�������仯 ����bid��binseq����ģ� ��ITEMID�ƺ�û�� �������˰�
        % �ص��Ǹ��������LU_VehType��BINSEQ��LU_VehType �⼸���ض��仯(PID/SID��Ҫ����) % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
        T{flaglastLUIdx,{'CoordLUBin','BINSEQ','LU_VehType'}} = T1{:,{'CoordLUBin','BINSEQ','LU_VehType'}};     chktLU(T);
        
    end

%% 3 ****************** ���ƽ��ѡ�� do2��˦βƽ�̣�/do3Array������ƽ�̣����� ��ȡ�޶��� output ******************
    if ISpingpu==1 && any(flagTiledArray~=0)  %���е�������ĳ��ƽ�̳ɹ���=1��=2��
        
        for ibin=1:length(do2Array) %do*Array ��������BIN
            
            if flagTiledArray(ibin)==0,continue; end  % ��ibinδƽ�� ����ѭ��
            
            if ISlastVehType==1 && flaggetSmallVeh == 1 && ibin == lastVehIdx, error('���һ���ȸ�������,����Ҫƽ��;');   end
            
            if flagTiledArray(ibin)~=0
                
                if flagTiledArray(ibin) == 1    % ��ibin����ƽ�̳ɹ�
                    dd = do3Array(ibin);
                elseif flagTiledArray(ibin) == 2      % ��ibin˦βƽ�̳ɹ�
                    dd = do2Array(ibin);
                end
                
                % ��ȡdd->LU��Table��ʽ
                T23 = getTableLU(dd);    chktLU(T23) ;
            end
            
            % ��ǰtable�������ڸ�bin���߼���            
            flagTileLUIdx = T.BINID ==ibin;  % V1" % flagTileLUIdx = T{:,'BINID'}==ibin; (���ߵȼۣ�
            
            %% CHECK �ܱ�T �� �ӱ� T23���ڲ�����
            % 1 �Ƿ��ibin�ڵ�LU�������ϸ����; �Ƿ�T��ѡ�е�flagTileLUIdx����LU���ڸ�ibin��Ҳ���ϸ������
            if ~issorted(sort(T23.BINSEQ),'strictascend') || ~issorted(sort(T.BINSEQ(flagTileLUIdx,:)'),'strictascend')
                T.BINSEQ(flagTileLUIdx,:)';
                T23.BINSEQ';
                sort(T23.BINSEQ)';
                sort(T23.BINID)';
                do2Array(ibin).LU.LU_Bin;
                sort(T.BINSEQ(flagTileLUIdx,:))';
                error('1');
            end
            % 2 �Ƿ��ܱ�T��BINID����ͳһbinid���Ƿ��ӱ�T23���Ƿ�ͬһbinid���Ƿ�ibin����T��LU��Ӧ��binid��
            if ~isscalar(unique(sort(T.BINID(flagTileLUIdx,:))))  || ~isscalar(unique(sort(T23.BINID)))  || ibin~= unique(sort(T.BINID(flagTileLUIdx,:)))
                sort(T23.BINID)'
                unique(sort(T23.BINID))
                unique(sort(T.BINID(flagTileLUIdx,:)))
                error('11');
            end
            % 3 �Ƿ�T23����������bin��LU����
            if sum(flagTileLUIdx) ~= height(T23)
                error('111');
            end
            % 4 �Ƿ��ܱ�T�ڳ�ʼOPID�뷵���ֱ���T23��OPID��ͬ;
            a = T.OPID(flagTileLUIdx,:)';
            b = T23.OPID';
            if  ~isequal(a,b)
                error('1');
            end
            
            %% ��Щ��仯?? ��ƽ��/gap�������bin�ڵ�LU�Ĳ�������(����,˳���)���ص��ܱ�T�� 5555
            % ĳ��bin�ڵ���,��binIDһ������仯; ��bin��LU_VehTypeһ������仯��
            % ��LID/Weight/LWHӦ�ò���仯; SID/PID��仯; ��ΪOPID OSID OID��ԭ�� idExchange����
            % ĳ��bin�ڵ���,��BINSEQ,CoordLUBinһ�������仯 ����bid��binseq����ģ�
            % ��LU_Itemһ����仯 ���ƺ�û���� �������˰�
            % �ص��Ǹ������� % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
            
            % ��ͼ����ƽ�̵��Ǹ�BIN��ͼ
            %     plotSolutionT(T23,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));
            
            % fixeme ����ʹ��ʱ����ϸ���� ��Щ�����仯�ĵ���
            T{flagTileLUIdx,{'CoordLUBin','BINSEQ','LU_Item'}} = ... %��������LU_Item�����л�,��Ȼ��;����,�����ᱨchktLU����.
                T23{:,{'CoordLUBin','BINSEQ','LU_Item'}};

            % 190421 ����
            T.LU_Bin = [T.BINID,T.BINSEQ]; % 190421 ����BINID/BINSEQ�ķ���LU_Bin,������ͼ; ����BINSEQҲ��
            % if parGap BINSEQ���ܱ仯 -> ���¼���BINSEQ ID����
            
            % �������Gap�ҳɹ��滻Gap,��bin�ڵ�LWH��RotaedҲҪ�滻��������T��
            if parGap % && flagGap Gap�����Ǳ���  18-12 �˴����ӶԻ�װ��϶�Ĵ��� - ����Table��ʽ
                T{flagTileLUIdx,{'LWH','Rotaed'}} = T23{:,{'LWH','Rotaed'}};
%                 T{flagTileLUIdx,{'BINSEQ'}} = T23{:,{'BINSEQ'}};
            end
            
            %        sortrows(T.LU_Item)'
            %        if ITEMID��ͬ ������CoordLUBin�ĳ��������ͬ ��ͬITEMID�����������Ա�
            %        ��table��ʽLU����chk ��Ҫ������
            %        chktLU(T) %�����޷�ͨ���Ŀ�����; ��LU_ItemӰ�첻��,������ע�� TODO
            
        end % END OF BIN
        
    end % END OF PINGPU

%% 4 ****************** ��Ի�ȡ��T ������󷵻ص�output���� ******************
    %���� ShowSEQ(������/��Ӧ�̺�/LID������ʾ����) �� tblorder -> ����ǿʹ��
    [T.ShowSEQ, T.tblorder] = getBBASeqTLU(T); 
    %       chktLU(T)  %      ���治ͨ��,������LU_Itemδ��ʱ����,�ں���˦βƽ�̺�. TODO

    % NOTE : �ڴ�֮ǰ��δ�ı�LU��˳�򣬸ı�LU˳��Ϊ������ʾ˳��ShowSEQ����
    [~,T.ttt] = sort(T.tblorder);
    T = sortrows(T,'ttt');                               % T = sortrows(T,'BINID')  % T.LID % T.BINID % T.BINSEQ
end







%% V1 V2ɾ�������ƺ�����������ʵ���õ�ע��
% % function T = HBinCombine(do,flaggetSmallVeh,do1,flagTiledArray,do2Array,do3Array)
% % global ISlastVehType ISpingpu parGap
% % %% 1 ******************��ȡչʾ˳�� do���� T=d.LU����ShowSEQ
% % T = getTableLU(do);
% % chktLU(T)
% %
% % % ��ͼ��ԭʼ��ƽ�̵�ͼ
% % %     plotSolutionT(T,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));
% %
% % %%
% % % ��1���������ڳ��ͺ�(���뻻)    ��2���������ڳ����(���,���ܻ�,���ʹ�) ��3�����̳��ڰ���˳��(���뻻) ��4������SID��Ӧ�̱��(�����,���ñ�??)
% % % ��5������ID�ͺ�LID(�����,���ñ�?) ��6�����̶Ѷ����ITEM(���,���ܻ�,���ʹ�,��;?) ��7�������㲿�����PID(�����,���ñ�?) ������8: չʾ˳��(���뻻)
% %
% % %% 2 ****************** ��Գ��ͱ仯 do1���� ��ȡ�޶��� output ******************
% % if ISlastVehType==1 && flaggetSmallVeh == 1 %���е������ҳ����滻�ɹ�
% %     T1 = getTableLU(do1);
% %     chktLU(T1)
% %
% %     % �滻T�е����һ���Ĳ������� ����T1
% %     lastVehIdx = max(T{:,'BINID'});
% %     flaglastLUIdx = T{:,'BINID'}==lastVehIdx;
% %
% %        %% ��Щ��仯??
% %        % ĳ��bin�ڵ���,��binIDһ������仯;
% %        % ��LID/Weight/LWHӦ�ò���仯; SID/PID��仯; ��ΪOPID OSID OID��ԭ��
% %        % ĳ��bin�ڵ���,��BINSEQ,CoordLUBin,LU_VehTypeһ�������仯 ����bid��binseq����ģ� ��ITEMID�ƺ�û�� �������˰�
% %        % �ص��Ǹ��������LU_VehType��BINSEQ��LU_VehType �⼸���ض��仯(PID/SID��Ҫ����) % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% %     T{flaglastLUIdx,{'CoordLUBin','BINSEQ','LU_VehType'}} = T1{:,{'CoordLUBin','BINSEQ','LU_VehType'}};
% %     chktLU(T)
% %     %%
% %
% %         % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% %                 % %     T{flaglastLUIdx,{'LU_VehType','BINSEQ','ShowSEQ'}} = ...
% %                 % %         T1{:,{'LU_VehType','BINSEQ','ShowSEQ'}};
% %                 % %     T{flaglastLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% %                 % %         T1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% %
% % %      [TLAST_Coord,TLAST_LWH,TLAST_Seq] = getReturnBBA(do1); %% ���з��ش���
% % %     TLAST_Seq1 = finalCheck( [TLAST_Coord,TLAST_LWH,TLAST_Seq],TLUIN_LAST); %����1�����㣻 ����2��ԭʼ.
% %
% % %     %����order�ı���,�˴��������һ��bin�����������޸�
% % %     lastVehIdx = max(T_Seq{:,'BINID'});
% % %     flaglastLUIdx = T_Seq1{:,'BINID'}==lastVehIdx;
% % %
% % %     % �ص��Ǹ�������ͳ����
% % %     T_Coord{flaglastLUIdx,:} = TLAST_Coord{:,:};
% % %     T_LWH{flaglastLUIdx,:} = TLAST_LWH{:,:};
% % %     % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% % %     T_Seq1{flaglastLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% % %         TLAST_Seq1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% %
% % %     T_Seq = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %����1�����㣻 ����2��ԭʼ.
% %
% % %     output_CoordLUBin(:,flaglastLUIdx) = TLAST_Coord;
% % %     output_LU_LWH(:,flaglastLUIdx) = TLAST_LWH;
% %
% %     % ������ز�����Ҫ�滻: ��1���������ڳ��ͺ� ��3�����̳��ڰ���˳�� ��8: չʾ˳��
% %         %     output_LU_Seq([1,3,4,5,7],flaglastLUIdx) = output_LU_Seq2([1,3,4,5,7],:); %[1,3,4,5,7]��ʾ���޸�����ļ���
% % %     output_LU_Seq([1,3,8],flaglastLUIdx) = TLAST_Seq([1,3,8],:); %[1,3,4,5,7,8]��ʾ���޸�����ļ���
% % %     output_LU_Seq([1,3,4,5,7,8],flaglastLUIdx) = TLAST_Seq([1,3,4,5,7,8],:); %[1,3,4,5,7,8]��ʾ���޸�����ļ���
% %
% %         %     i = [4,5,7]; % ��ЩӦ���ǲ�����
% %         %     if sum(output_LU_Seq(i,flaglastLUIdx) ~= output_LU_Seq2(i,:) ) >0, error('�����ı���, ����'); end
% % end
% %
% % %% 3 ****************** ���ƽ��ѡ�� do2��˦βƽ�̣�/do3Array������ƽ�̣����� ��ȡ�޶��� output ******************
% % if ISpingpu==1 && any(flagTiledArray~=0)  %���е�������ĳ��ƽ�̳ɹ���=1��=2��
% %
% %     for ibin=1:length(do2Array) %do*Array ��������BIN
% %         if flagTiledArray(ibin)==0  % ��ibinδƽ�� ����ѭ��
% %             continue;
% %         end
% %         if ISlastVehType==1 && flaggetSmallVeh == 1 && ibin == lastVehIdx
% %             error('���һ���ȸ�������,����Ҫƽ��;');
% %         end
% %
% %         if flagTiledArray(ibin)==1    % ��ibin����ƽ�̳ɹ�
% %             dd = do3Array(ibin);
% %
% %         % ������װ��϶�Ż�     %             if parGap==1,   [dd.LU] = getMixedGap(dd.LU, dd.Veh);    end
% %
% %             % ��ȡLU��Table��ʽ
% %             T23 = getTableLU(dd);   %T23 = getTableLU(do3Array(ibin));
% %             chktLU(T23) ;
% %         end
% %         if flagTiledArray(ibin)==2      % ��ibin˦βƽ�̳ɹ�
% %             dd = do2Array(ibin);
% %
% %         % ������װ��϶�Ż�     %             if parGap==1,   [dd.LU] = getMixedGap(dd.LU, dd.Veh);    end
% %
% %             % ��ȡLU��Table��ʽ
% %             T23 = getTableLU(dd);     % T23 = getTableLU(do2Array(ibin));
% %             chktLU(T23);
% %         end
% %
% %         % ��ǰtable�������ڸ�bin���߼���
% %         flagTileLUIdx = T{:,'BINID'}==ibin;
% %
% %        %% CHECK 1 �Ƿ��ibin�ڵ�LU�������ϸ����; �Ƿ�T��ѡ�е�flagTileLUIdx����LU���ڸ�ibin��Ҳ���ϸ������
% %        if ~issorted(sort(T23.BINSEQ),'strictascend') || ~issorted(sort(T.BINSEQ(flagTileLUIdx,:)'),'strictascend')
% %                T.BINSEQ(flagTileLUIdx,:)'
% %                T23.BINSEQ'
% %                sort(T23.BINSEQ)'
% %                sort(T23.BINID)'
% %                do2Array(ibin).LU.LU_Bin
% %                sort(T.BINSEQ(flagTileLUIdx,:))'
% %                error('1');
% %        end
% %        % 2 �Ƿ��ܱ�T��BINID����ͳһbinid���Ƿ��ӱ�T23���Ƿ�ͬһbinid���Ƿ�ibin����T��LU��Ӧ��binid��
% %        if ~isscalar(unique(sort(T.BINID(flagTileLUIdx,:))))  || ~isscalar(unique(sort(T23.BINID)))  || ibin~= unique(sort(T.BINID(flagTileLUIdx,:)))
% %             sort(T23.BINID)'
% %             unique(sort(T23.BINID))
% %             unique(sort(T.BINID(flagTileLUIdx,:)))
% %             error('11');
% %        end
% %        % 3 �Ƿ�T23����������bin��LU����
% %        if sum(flagTileLUIdx) ~= height(T23)
% %            error('111');
% %        end
% %        % 4 �Ƿ��ܱ�T�ڳ�ʼOPID�뷵���ֱ���T23��OPID��ͬ;
% %        a = T.OPID(flagTileLUIdx,:)';
% %        b = T23.OPID';
% %        if  ~isequal(a,b)
% %            error('1');
% %        end
% %
% %        %% ��Щ��仯??
% %        % ĳ��bin�ڵ���,��binIDһ������仯; ��bin��LU_VehTypeһ������仯��
% %        % ��LID/Weight/LWHӦ�ò���仯; SID/PID��仯; ��ΪOPID OSID OID��ԭ�� idExchange����
% %        % ĳ��bin�ڵ���,��BINSEQ,CoordLUBinһ�������仯 ����bid��binseq����ģ�
% %        % ��LU_Itemһ����仯 ���ƺ�û���� �������˰�
% %        % �ص��Ǹ������� % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% % %        sortrows(T.LU_Item)'
% % %        sortrows(T23.LU_Item)'
% % %        T.LU_Item
% % %        T{flagTileLUIdx,{'CoordLUBin','BINSEQ'}} = T23{:,{'CoordLUBin','BINSEQ'}};
% %
% %     % ��ͼ����ƽ�̵��Ǹ�BIN��ͼ
% % %     plotSolutionT(T23,struct2table(structfun(@(x) x',d.Veh,'UniformOutput',false)));
% %
% %         T{flagTileLUIdx,{'CoordLUBin','BINSEQ','LU_Item'}} = ... %��������LU_Item�����л�,��Ȼ��;����,�����ᱨchktLU����.
% %             T23{:,{'CoordLUBin','BINSEQ','LU_Item'}};
% %
% %         % �������Gap�ҳɹ��滻Gap,��bin�ڵ�LWH��RotaedҲҪ�滻��������T��
% %         if parGap % && flagGap Gap�����Ǳ���
% %             T{flagTileLUIdx,{'LWH','Rotaed'}} = ... %��������LU_Item�����л�,��Ȼ��;����,�����ᱨchktLU����.
% %                 T23{:,{'LWH','Rotaed'}};
% %         end
% %
% % %        sortrows(T.LU_Item)'
% % %        chktLU(T) %�����޷�ͨ���Ŀ�����; ��LU_ItemӰ�첻��,������ע�� TODO
% %
% %             %% �����Ǵ��
% %             %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','ShowSEQ'}} = ...
% % %            T23{:,{'LU_VehType','BINSEQ','ShowSEQ'}};
% % %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% % %            T23{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% % %        T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH'}} = ...
% % %            T23{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH'}};
% %
% %
% %             %     T{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH','BINID','ITEMID','Weight'}} = ...
% %             %         T2{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ','CoordLUBin','LWH','BINID','ITEMID','Weight'}};
% %
% %             % output_LU_Seq=T{:,{'','BINID',',''ITEMID','Weight'}}'
% %
% %             % % %     [TPP_Coord,TPP_LWH,TPP_Seq]= getReturnBBA(do2Array(ibin)); %% ���з��ش���
% %             % % %     %     [output_CoordLUBin2,output_LU_LWH2,output_LU_Seq2]= getReturnBBA1(do2Array(ibin)); %% ���з��ش���
% %             % % %     TPP_Seq1 = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP1);
% %             % % % %     if  sum(flagTiled(:)~=0)>1
% %             % % % %         sum(flagTiled(:)~=0)
% %             % % % %         1
% %             % % % %     end
% %             % % % %     if flagTiled(ibin)==1
% %             % % % %         if sum(flagTiled(:)==1)>1
% %             % % % %             sum(flagTiled(:)==1)
% %             % % % %             1
% %             % % % %         end
% %             % % % %         TPP_Seq = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP1); %����1�����㣻 ����2��ԭʼ.
% %             % % % %     end
% %             % % % %     if flagTiled(ibin)==2
% %             % % % %         if  sum(flagTiled(:)==2)>1
% %             % % % %             sum(flagTiled(:)==2)
% %             % % % %             1
% %             % % % %         end
% %             % % % %         TPP_Seq = finalCheck([TPP_Coord,TPP_LWH,TPP_Seq],TLUIN_PP2); %����1�����㣻 ����2��ԭʼ.
% %             % % % %     end
% %             % % %
% %             % % %     % �ҳ�ƽ��ibin�����е������߼�ֵ
% %             % % %     flagTileLUIdx = T_Seq1{:,'BINID'} == ibin;
% %             % % %
% %             % % %     % �ص��Ǹ�������ͳ����
% %             % % %     T_Coord{flagTileLUIdx,:} = TPP_Coord{:,:};
% %             % % %     T_LWH{flagTileLUIdx,:} = TPP_LWH{:,:};
% %             % % %     % LU_VehType   'BINID'   BINSEQ   SID    LID    'ITEMID'    PID  ShowSEQ   'Weight'
% %             % % %     T_Seq1{flagTileLUIdx,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}} = ...
% %             % % %         TPP_Seq1{:,{'LU_VehType','BINSEQ','SID','LID','PID','ShowSEQ'}};
% %
% %
% %             %         T_Seq = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %����1�����㣻 ����2��ԭʼ.
% %
% %             % ������ز�����Ҫ�滻: ��1���������ڳ��ͺ� ��3�����̳��ڰ���˳�� ��8: չʾ˳��
% %             %     output_LU_Seq([1,3,8],flagTileLUIdx) = output_LU_Seq2([1,3,8],:);
% %             %     output_LU_Seq([1,3,4,5,7,8],flagTileLUIdx) = output_LU_Seq2([1,3,4,5,7,8],:); %[1,3,4,5,7,8]��ʾ���޸�����ļ���
% %
% %             %     i = [4,5,7]; % ��ЩӦ���ǲ�����
% %             %     if sum(output_LU_Seq(i,flagTileLUIdx) ~= output_LU_Seq3(i,:) ) >0,
% %             %          output_LU_Seq(i,flagTileLUIdx) ~= output_LU_Seq3(i,:) ;         output_LU_Seq(i,flagTileLUIdx);         output_LU_Seq3(i,:) ;
% %             %         warning('�����ı���, ����'); end
% %     end
% % end
% %
% % %% 18-12 �˴����ӶԻ�װ��϶�Ĵ��� - ����Table��ʽ
% %
% % % ****************** ��Գ���ѡ�� ��ȡ�޶��� output ******************
% %
% % T = getShowSeq(T); %����ShowSEQ tblorder
% %
% % %        chktLU(T)  %      ���治ͨ��,������LU_Itemδ��ʱ����,�ں���˦βƽ�̺�. TODO
% %
% %
% % %% tableתΪ�ṹ����ж��Ƿ���������
% % % lu = table2struct(T,'ToScalar',true)
% % % lu = (structfun(@(x) x',lu,'UniformOutput',false));
% %
% %     %% ����T������LU_LUinBin���������ж�
% % %     T2 = sortrows(T,{'BINID','BINSEQ'},{'ascend','ascend'})
% %     %     x=T(T.BINID==2&T.ID==1,{'ID','LID','PID','H','Weight','CoordLUBin','BINSEQ','ShowSEQ','ITEMID','ITEMSEQ'})
% % %     x=T2(:,{'ID','LID','PID','H','Weight','X','Y','Z','ITEMID','ITEMSEQ','BINID','BINSEQ','ShowSEQ'})
% %     % if ITEMID��ͬ ������CoordLUBin�ĳ��������ͬ ��ͬITEMID�����������Ա�
% %     % ��table��ʽLU����chk ��Ҫ������
% %
% %
% %
% % % NOTE : �ڴ�֮ǰ��δ�ı�LU��˳��
% % % ����BBA�����ʽ��JAR
% % [~,T.ttt] = sort(T.tblorder);
% % T = sortrows(T,'ttt');
% % % T = sortrows(T,'BINID')
% % % T.LID
% % % T.BINID
% % % T.BINSEQ
% % end
% %
% %
% % %% ע�ͺ�
% %
% %
% % %% USELESS
% % % [T_Coord,T_LWH,T_Seq] = getReturnBBA(daBest(bestOne)); %���ж��,���ص�һ�����Ž�
% % % T_Seq.tblorder
% % % T_Seq1 = T_Seq;
% %
% % % T.SID = oD.SID;
% % % T.PID = oD.PID;
% % % T.LID = oD.LID;
% % % Tseq = T(:,{'LU_VehType','BINID','BINSEQ','SID','LID','ITEMID','PID','ShowSEQ','Weight','tblorder'});
% %
% % % T_Seq1 = finalCheck([T_Coord,T_LWH,T_Seq],TLUIN); %����1�����㣻 ����2��ԭʼ.
% %
% % %  [output_CoordLUBin,output_LU_LWH,output_LU_Seq]= getReturnBBA1(daBest(bestOne)); %% ���з��ش���
