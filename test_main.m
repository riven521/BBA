% 测试函数:调用主函数
clc;  clear;  close all; 

%     LUID = [5236934585 4 5236934585 5236934585];
%     BINLWH = [8;10;4];  %[6; 10; 4];
%     LULWH = [3 2 3 3; 6 4 6 6; 2 2 1 2];    
%     PARA = 0;    
%     LUBUFF = [0,0];
%     BINBUFF = [0,0,0];       
    
% 企业测试版1
    LUID = [2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2] ;
            LUSID = ones(size(LUID));
        LUPID = ones(size(LUID));
        LUISROTA = ones(size(LUID));
        LUWEIGHT = ones(size(LUID));
        LUMARGIN = zeros(4,size(LUID,2));
         VEHWEIGHT = [1000];
    BINLWH = [2380, 9500, 2600];  %[6; 10; 4];
    LULWH = [1955 760 550 ;
1220 820 1100 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1955 760 550 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1220 820 1100 ;
1955 760 550 ]';    

%     PARA = [ 1 0] ;    
%     LUBUFF = [100;100];
%     BINBUFF = [0;0;0];       
    BINID = [1];
% %     LUID = [1 2 1 2];
% %     BINLWH = [10;20;4]';
% %     LULWH = [2 5 4;  7 5 4; 2 5 2; 7 5 4];
% %     LUBUFF = [0,0];
% %     BINBUFF = [0,0,0];
    
    
    % whichStripH 1 best 2 first 3 next; 
    % whichBinH 1 best; TODO 增加其它分批方式
    % whichSortItemOrder 1 长高递减 2 最短边递减; 
    % whichRotation 1:允许rotation 0:禁止 
    % whichRotationHori 0-2
    % whichRotationAll 0 1：全部以Item rotation后的为主
    % whichRotationBin 0 1：全部以Bin rotation后的为主
    % 0 1 2
    
    %% 不同参数多次计算
%                             [Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq] = ...
%                             BBA_Main(LUID,LULWH,BINID,BINLWH,BINBUFF);
                                                    [Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq] = ...
                            BBA_Main(LUID,LULWH,BINID,BINLWH,...
                            LUSID, LUPID, LUISROTA, LUMARGIN, LUWEIGHT, VEHWEIGHT  );
% %      n=5;d=getRandDa(n);save('rndDa.mat','d');
% %      load('rndDa.mat')
% %          LUID =d.LUArray.ID ; 
% %     LULWH=d.LUArray.LWH ;   %LWHREAL 真实尺寸    
% %     BINLWH =d.BinArray.LWH ; %LWHREAL 真实尺寸    
% %     % 增加间隙 -
% %     LUBUFF= d.LUArray.BUFF ; %BUFF 托盘LU的间隙
% %     BINBUFF=d.BinArray.BUFF ; %BUFF 车辆BIN的间隙

     
% % % res = zeros(1,36);
% % % respara = zeros(7,36);
% % % r = 1;
% % % for i = 1:3 %1-3 best first next均可
% % %     for j=1:2 %1-2 排序可多选
% % %         for k=1:1 %0-1 默认1可旋转
% % %             for l=1:2 %0-2 0不好
% % %                 for p =0:1
% % %                     for o = 0:0 %车辆rota不如物流rota
% % %                         PARA = [i 1 j k l p o];
% % %                         figure('name',num2str(PARA));
% % %                         [Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq,d] = ...
% % %                             BBA_Main(LUID,LULWH,BINLWH,PARA,LUBUFF,BINBUFF);
% % %                         res(1,r) = mean(d.Strip.ItemloadingrateLimit); %Itemloadingrate
% % %                         respara(:,r) = PARA';
% % %                         r=r+1;
% % %                     end
% % %                 end
% % %             end
% % %         end
% % %     end
% % % end
% % % 
% % % close all
% % % [x,y]=max(res);
% % % [a,~]=find(res(:)==x);
% % % for b=1:length(a)
% % %     PARA=respara(:,a(b));
% % %     figure('name',num2str(PARA)); %figure从1开始i排序
% % %     [Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq,d] = ...
% % %         BBA_Main(LUID,LULWH,BINLWH,PARA,LUBUFF,BINBUFF);
% % % end


% % PARA
% % y
% % x
% % res