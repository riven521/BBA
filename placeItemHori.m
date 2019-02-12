function [flag]  =placeItemHori(LWH,isRota,margin,HoriOrVert)
% placeItemHori ==> ����LU��Rotaed���( LU��Ż����� ȴ���ھ���߽�Ŀ�� �ĸ�Сѡ�ĸ� 
%   ��ÿ�ν���RunAlgorithmǰ���� -> ��ȡLU�ĵ�ǰ״̬

        % ��ȡ����margin��LWH
        [LWHbuff] = LWHbuffer(LWH, margin);
        
        % ��ȡ�Ƿ���ת�ı��flag
        flag1 = LWHbuff(1,:) < LWHbuff(2,:);  %��<��/��,       1����Horizontal 0����Horizontal
        flag2 = isRota == 1;                                    % 1: ������ת 0 ������
        
        if HoriOrVert == 1 %��Horizontal �ڷ�
            error('111111111111');
            flag = flag1 & flag2;  %������ת��Ŀǰ��Horizontal��ʽ�ڷ�
        elseif HoriOrVert == 0 %��Vertical �ڷ�
            error('111111111111');
            flag = ~flag1 & flag2;  %������ת��Ŀǰ��Horizontal��ʽ�ڷ�
        elseif HoriOrVert ==2
            error('111111111111');
            flag = false(size(flag1)); %��������ԭ�ⲻ��
        else
            
            wVeh = HoriOrVert;
            x=mod(wVeh,LWHbuff(1,:));   %��������ȶ�Hori�ŷ�ȡ����
            y=mod(wVeh,LWHbuff(2,:));   %��������ȶ�Vert�ŷ�ȡ����
            
            flag3 = x <= y; %�Ƿ�Hori�ŷű�Vert�ŷ�������С 1 ϣ��Hori�ڷ� 0 ����
            fxor = ~xor(flag1,flag3); %XOR ���Ա�           
            flag =  fxor & flag2;
            
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