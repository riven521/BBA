function [] = plotSolutionBBA(c,lwh,m,d)
% plotSolutionBBA(output_CoordLUBin,output_LU_LWH,output_LU_Seq,daBest(bestOne).Veh);
% ��ͼ����:��άBPP
Veh = d.Veh;
LU = d.LU;
Item = d.Item;
Type = 'LU';
% plot3DStrip(LU,Item,Veh,Type);
% % չʾ�ò��� m: 1-9: [LU_Bin(1,:); LU_Bin(2,:); SID; LID; PID; LU_Item(1,:); LU_Item(2,:); hLU; LU_VehType];
% m : tmpShow =[9,1,2,3,4,6,5];  + չʾ˳��  
% m1 LU_VehType  m2 LU_Bin(1,:) m3 LU_Bin(2,:) m4  SID; m5 LID;  m6
% LU_Item(1,:); m7 PID;  m8 չʾ˦β m9 չʾ˳��

global ISplotPause ISplotShowType
if ISplotShowType == 1
    % 1: ����LID������ɫ
    nIDType = unique(m(5,:));
elseif ISplotShowType == 2
    % 2: ����˦β���������ɫ
    nIDType = unique(m(8,:));
end

nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
nBin = max(m(2,:)); %bin�ĸ���;

figure('name',num2str([nBin]));
j = 1;
tmp = m(8,1);
for ibin=1:nBin

    subplot(2,ceil((nBin)/2),ibin);
    
    f = m(2,:)==ibin;
    yxz = lwh(:,f);
    coord = c(:,f);
    
    if ISplotShowType == 1
   % 1: ����LID������ɫ 
   lid  = m(5,f);
   elseif ISplotShowType == 2
   % 2: ����˦β���������ɫ
   lid  = m(8,f);
    end

    % LU��Bin�ڵ�˳��  LU_Bin(2,:)
    seq = m(3,f);

    for i=1:numel(seq)
        idx =find(seq==i); % LU����Bin�ڵĻ�ͼ˳��
        LUColor = 0.8*nColors(nIDType==lid(idx), : );

        % �����̲���������ʱ, �ж�
        if j~=1 && tmp ~=m(8,j) 
                 if ISplotPause>=0 ,      pause(ISplotPause);   end
        end        
        tmp = m(8,j);
        j=j+1;
        
        % plot
        plotcube(yxz(:,idx)',coord(:,idx)',0.7,LUColor);
        
        % Set the lable and the font size
        axis equal;         grid on;        view(111,33); %view(60,40);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        
        % �����ĳ���ߵ��������ʵĳ���
        xlim([0 Veh.LWH(1,unique(m(1,f)))]);         ylim([0 Veh.LWH(2,unique(m(1,f)))]);       zlim([0 Veh.LWH(3,unique(m(1,f)))]);        
    end
end

end
