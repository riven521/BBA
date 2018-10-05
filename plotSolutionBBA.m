function [] = plotSolutionBBA(c,lwh,m,Veh)
% ��ͼ����:��άBPP

% 1: ����LID������ɫ
 nIDType = unique(m(5,:));

nColors = hsv(length(nIDType)); %��ͬ����LU���費ͬ��ɫ
nBin = max(m(2,:)); %bin�ĸ���;

figure('name',num2str([nBin]));
for ibin=1:nBin

    subplot(2,ceil((nBin)/2),ibin);
    
    f = m(2,:)==ibin;
    yxz = lwh(:,f);
    coord = c(:,f);
    
   lid  = m(5,f);

    % LU��Bin�ڵ�˳��
    seq = m(3,f);
    
    for i=1:numel(seq)
        idx =find(seq==i); % LU����Bin�ڵĻ�ͼ˳��
        LUColor = 0.8*nColors(nIDType==lid(idx), : );
        plotcube(yxz(:,idx)',coord(:,idx)',0.7,LUColor);
        
        % Set the lable and the font size
        axis equal;         grid on;        view(111,33); %view(60,40);
        xlabel('X','FontSize',10);         ylabel('Y','FontSize',10);         zlabel('Z','FontSize',10);
        
        % �����ĳ���ߵ��������ʵĳ���
        xlim([0 Veh.LWH(1,unique(m(1,f)))]);         ylim([0 Veh.LWH(2,unique(m(1,f)))]);       zlim([0 Veh.LWH(3,unique(m(1,f)))]);        
    end
end

end
