function [flaggetSmallVeh,do1] = HBinChange(maind,do,p)
% HBinChange ==> 更换车型

    % 2.1 初始化
    allidxVehType = length(unique(do.Veh.ID)); %此算例车型数量(未排除相同车型)
    flaggetSmallVeh = 0;
    
                                                                        % 准备替换d1 %  对调Lu.LWH的长宽 -< 之前是宽长 (已放入getdinLastVeh中)
                                                                        %     d1.LU.LWH([1,2],:) = flipud(d1.LU.LWH([1,2],:)); 
                                                                        %d1 = getdinLastVeh(do);
    % 2.2 CHANGE 从maind中拿输入数据
    luIdx = do.LU.LU_Bin(1,:) == max(do.LU.LU_Bin(1,:));
    d1 = getPartdinThisBin(maind,luIdx);                    
                    
    % 2.3 当车辆类型多于2种,才进行替换
    while(allidxVehType>=2)
        % 2.1 获取最后车型并运行算法 % 从最后一辆车不断往前循环; until第二辆车; 此处假设
        d1.Veh = structfun(@(x) x(:,allidxVehType), do.Veh,'UniformOutput',false); %从最后一种车型开始考虑
                                            %disp(d1.Veh.LWH)
                                            %TLUIN_LAST = struct2table(structfun(@(x) x',d1.LU,'UniformOutput',false));
        
        [d1.LU] = setLULWHwithbuff(d1.LU, d1.Veh);                                                                                                                   %d1 = RunAlgorithm(d1,pA(iAlg));
        do1 = RunAlgorithm(d1,p);   %针对少数的最后一个Bin的输入lastd进行运算 555555555555555555555
                                            %     plotSolution(do1,pA(iAlg));

                            %             do1.LU.LU_VehType = ones(size(do1.LU.ID)) * do1.Veh.order(1); % 针对车型选择,增加变量LU_VehType : 由于Veh内部按体积递减排序,获取order的第一个作为最大值
        % do1.LU.LU_VehType = ones(size(do1.LU.ID))*do.Veh.order(allidxVehType); % 补充变量LU_VehType
        [do1.LU,do1.Item] = setLCwithoutbuff(do1.LU,do1.Item);

        % 2.4 判断该车型是否可用
        % 由于Veh内部按体积递减排序,获取order的第个作为当前对应真车型索引号
        % 判断: 是否改为第allidxVehType(小)车型后,1个车辆可以放下;
        if max(do1.LU.LU_Bin(1,:)) == 1
            fprintf(1,'       Exsiting 车型更换 in HBinChange (do1)...\n');
                                                                                                                               %do1.LU.LU_VehType = ones(size(do1.LU.ID))*do.Veh.order(allidxVehType); % 补充变量LU_VehType
            flaggetSmallVeh=1;
            break;
        end
        
        % 2.5 若放不下,选择更大车型 -> allidxVehType递减 do1.Veh赋予空值
        allidxVehType= allidxVehType-1;
    end
    
    %% chk 2 运行车型调整算法,
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