% ��Item����ˮƽ����ֱ����: LU��Item������ʹ��
function [flag]  =placeItemHori(LWH,isRota,HoriOrVert)
        % �����ж��Ƿ���ITEM�ı仯, �����, Ҫ���Ӷ�LU�仯��Ӱ��
        
        % ��ȡ�Ƿ���ת�ı��flag
        flag1 = LWH(1,:) < LWH(2,:);  %��<��/��, 1����Horizontal 0����Horizontal
        flag2 = isRota == 1;                       % 1: ������ת 0 ������
        if HoriOrVert == 1 %��Horizontal �ڷ�
            flag = flag1 & flag2;  %������ת��Ŀǰ��Horizontal��ʽ�ڷ�
        elseif HoriOrVert == 0 %��Vertical �ڷ�
            flag = ~flag1 & flag2;  %������ת��Ŀǰ��Horizontal��ʽ�ڷ�
        else
            flag = false(size(flag1)); %��������ԭ�ⲻ��
        end
end
    

% % % ��Item����ˮƽ����ֱ����: LU��Item������ʹ��
% % function [L, flag]  =placeItemHori(Item,HoriOrVert)
% %         % ��ȡ�Ƿ���ת�ı��flag
% %         flag1 = Item.LWH(1,:) < Item.LWH(2,:);  %��<��/��, 1����Horizontal 0����Horizontal
% %         flag2 = Item.isRota == 1;                       % 1: ������ת 0 ������
% %         if HoriOrVert == 1 %��Horizontal �ڷ�
% %             flag = flag1 & flag2;  %������ת��Ŀǰ��Horizontal��ʽ�ڷ�
% %         elseif HoriOrVert == 0 %��Vertical �ڷ�
% %             flag = ~flag1 & flag2;  %������ת��Ŀǰ��Horizontal��ʽ�ڷ�
% %         else
% %             flag = false(size(flag1)); %��������ԭ�ⲻ��
% %         end
% % 
% %         %�ȸ�ֵ�������Ҫ��Item��������            
% %         L            = Item.LWH;
% %         L(1,flag) = Item.LWH(2,flag);
% %         L(2,flag) = Item.LWH(1,flag);
% % end
% %    

% % function [ItemHori,flag] = horiOrient(Item,isRota)
% % % ��ͨ����; ��Itemת��ΪHorizontal Orientation�ڷţ�����ƽ����ײ���
% % % NOTE : ��ǰ��ֻ�ı���Ҫ�ı�ģ�Ŀǰ��ȫ���ı� ItemHori = zeros(size(Item));    ItemHori = Item;    
% %  
% %     ItemHori = zeros(size(Item));
% %     ItemHori(1,:) = max(Item);
% %     ItemHori(2,:) = min(Item);
% %     
% %     flag = Item(1,:) < Item(2,:);
% %     flag = double(flag);             %logicalת��Ϊdouble
% %  
% % %     tmpminv = min(Item);
% % %     tmpmaxv = max(Item);
% % %     ItemHori(1,flag) = tmpmaxv(1,flag);   
% % %     ItemHori(2,flag) = tmpminv(1,flag);    
% % end