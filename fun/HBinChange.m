function [flaggetSmallVeh,do1] = HBinChange(maind,do,p)
% HBinChange ==> ��������

    % 2.1 ��ʼ��
    allidxVehType = length(unique(do.Veh.ID)); %��������������(δ�ų���ͬ����)
    flaggetSmallVeh = 0;
    
                                                                        % ׼���滻d1 %  �Ե�Lu.LWH�ĳ��� -< ֮ǰ�ǿ� (�ѷ���getdinLastVeh��)
                                                                        %     d1.LU.LWH([1,2],:) = flipud(d1.LU.LWH([1,2],:)); 
                                                                        %d1 = getdinLastVeh(do);
    % 2.2 CHANGE ��maind������������
    luIdx = do.LU.LU_Bin(1,:) == max(do.LU.LU_Bin(1,:));
    d1 = getPartdinThisBin(maind,luIdx);                    
                    
    % 2.3 ���������Ͷ���2��,�Ž����滻
    while(allidxVehType>=2)
        % 2.1 ��ȡ����Ͳ������㷨 % �����һ����������ǰѭ��; until�ڶ�����; �˴�����
        d1.Veh = structfun(@(x) x(:,allidxVehType), do.Veh,'UniformOutput',false); %�����һ�ֳ��Ϳ�ʼ����
                                            %disp(d1.Veh.LWH)
                                            %TLUIN_LAST = struct2table(structfun(@(x) x',d1.LU,'UniformOutput',false));
        
        [d1.LU] = setLULWHwithbuff(d1.LU, d1.Veh);                                                                                                                   %d1 = RunAlgorithm(d1,pA(iAlg));
        do1 = RunAlgorithm(d1,p);   %������������һ��Bin������lastd�������� 555555555555555555555
                                            %     plotSolution(do1,pA(iAlg));

                            %             do1.LU.LU_VehType = ones(size(do1.LU.ID)) * do1.Veh.order(1); % ��Գ���ѡ��,���ӱ���LU_VehType : ����Veh�ڲ�������ݼ�����,��ȡorder�ĵ�һ����Ϊ���ֵ
        % do1.LU.LU_VehType = ones(size(do1.LU.ID))*do.Veh.order(allidxVehType); % �������LU_VehType
        [do1.LU,do1.Item] = setLCwithoutbuff(do1.LU,do1.Item);

        % 2.4 �жϸó����Ƿ����
        % ����Veh�ڲ�������ݼ�����,��ȡorder�ĵڸ���Ϊ��ǰ��Ӧ�泵��������
        % �ж�: �Ƿ��Ϊ��allidxVehType(С)���ͺ�,1���������Է���;
        if max(do1.LU.LU_Bin(1,:)) == 1
            fprintf(1,'       Exsiting ���͸��� in HBinChange (do1)...\n');
                                                                                                                               %do1.LU.LU_VehType = ones(size(do1.LU.ID))*do.Veh.order(allidxVehType); % �������LU_VehType
            flaggetSmallVeh=1;
            break;
        end
        
        % 2.5 ���Ų���,ѡ������� -> allidxVehType�ݼ� do1.Veh�����ֵ
        allidxVehType= allidxVehType-1;
    end
    
    %% chk 2 ���г��͵����㷨,
    if flaggetSmallVeh==1
        % d1;
        % do1
%         t1 = struct2table(structfun(@(x) x',d1.LU,'UniformOutput',false));
%         to1 = struct2table(structfun(@(x) x',do1.LU,'UniformOutput',false));
        chkLUnewold(d1.LU,do1.LU);
        chktLU(do1.LU);
    else
        do1 = false;
    end

end