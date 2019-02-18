function [LWHBuff,OLWH] = LWHbuffer(LWH, margin, rota)
% add LWH��buffer ��������ת״̬
% ����3�����У�Ĭ��Ϊ�򲻻ָ���ת
    narginchk(2,3);

    OLWH = LWH;
    
    LWHBuff = LWH;    
    if nargin > 2
        rota = logical(rota);
        % ����ת�Ĳ���������buff������ת�ĳ����滻������buff��
        LWHBuff(1, rota) = LWH(2, rota) + (margin(1, rota) +margin(2, rota)); %��ȣ����ң�
        LWHBuff(1, ~rota) = LWH(1, ~rota) + (margin(1, ~rota) + margin(2, ~rota));
        
        LWHBuff(2, rota) = LWH(1, rota) + (margin(3, rota) +margin(4, rota)); %���ȣ����£�
        LWHBuff(2, ~rota) = LWH(2, ~rota) + (margin(3, ~rota) + margin(4, ~rota));
    else
        LWHBuff(1, :) = LWH(1, :) + (margin(1, :) +margin(2, :));  %��ȣ����ң�
        LWHBuff(2, :) = LWH(2, :) + (margin(3, :) + margin(4, :));  %���ȣ����£�
    end

end