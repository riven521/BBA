% ���Ժ���:����������
clc;  clear;  close all; 

%     LUID = [5236934585 4 5236934585 5236934585];
%     BINLWH = [8;10;4];  %[6; 10; 4];
%     LULWH = [3 2 3 3; 6 4 6 6; 2 2 1 2];    
%     PARA = 0;    
%     LUBUFF = [0,0];
%     BINBUFF = [0,0,0];       
    
% ��ҵ���԰�
% % %     LUID = [2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2] ;
% % %     BINLWH = [2380, 9500, 2600];  %[6; 10; 4];
% % %     LULWH = [1955 760 550 ;
% % % 1220 820 1100 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1955 760 550 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1220 820 1100 ;
% % % 1955 760 550 ];    
% % %     PARA = [0] ;    
% % %     LUBUFF = [0,0];
% % %     BINBUFF = [0,0,0];       

% %     LUID = [1 2 1 2];
% %     BINLWH = [10;20;4]';
% %     LULWH = [2 5 4;  7 5 4; 2 5 2; 7 5 4];
% %     LUBUFF = [0,0];
% %     BINBUFF = [0,0,0];
    
    
    % whichStripH 1 best 2 first 3 next; 
    % whichBinH 1 best; TODO ��������������ʽ
    % whichSortItemOrder 1 ���ߵݼ� 2 ��̱ߵݼ�; 
    % whichRotation 1:����rotation 0:��ֹ 
    % whichRotationHori 0-2
    % whichRotationAll 0 1��ȫ����Item rotation���Ϊ��
    % whichRotationBin 0 1��ȫ����Bin rotation���Ϊ��
    % 0 1 2
    
    %% ��ͬ������μ���
    
     n=5;da=getRandDa(n);save('rndDa.mat','da');
     load('rndDa.mat')
         LUID =da.LUArray.ID ; 
    LULWH=da.LUArray.LWH ;   %LWHREAL ��ʵ�ߴ�    
    BINLWH =da.BinArray.LWH ; %LWHREAL ��ʵ�ߴ�    
    % ���Ӽ�϶ -
  LUBUFF= da.LUArray.BUFF ; %BUFF ����LU�ļ�϶
    BINBUFF=da.BinArray.BUFF ; %BUFF ����BIN�ļ�϶

     
    res = zeros(1,36);
    respara = zeros(7,36);
    r = 1;
    for i = 1:3 %1-3 best first next����
        for j=1:2 %1-2 ����ɶ�ѡ
            for k=1:1 %0-1 Ĭ��1����ת
                 for l=0:2 %0-2 0����
                    for p =0:1
                        for o = 0:0 %����rota��������rota
                    PARA = [i 1 j k l p o];
                    figure('name',num2str(PARA));
                    [Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq,da] = ...
                        BBA_Main(LUID,LULWH,BINLWH,PARA,LUBUFF,BINBUFF);
                    res(1,r) = mean(da.StripArray.ItemloadingrateLimit); %Itemloadingrate
                    respara(:,r) = PARA';
                    r=r+1;
                        end
                    end
                end
            end
        end
    end
    close all    
    [x,y]=max(res);
    [a,~]=find(res(:)==x);
    for b=1:length(a)
        PARA=respara(:,a(b));
        figure('name',num2str(PARA)); %figure��1��ʼi����
        [Res1_LUBeBinMatrix,Res2_CoordLUBin,Res3_LWHRota,Res4_DrawSeq,da] = ...
            BBA_Main(LUID,LULWH,BINLWH,PARA,LUBUFF,BINBUFF);
    end
PARA
y
x
res