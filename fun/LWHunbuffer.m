function OLWH = LWHunbuffer(LWH, margin, rota)
% remove LWH��buffer �һָ���ת״̬
% ��û�в���3���򲻻ָ���ת

    narginchk(2,4);
    
    OLWH = LWH;     
    
    if nargin > 2
        rota = logical(rota); 
        LWH(1, rota ) = OLWH(2, rota);
        LWH(2, rota ) = OLWH(1, rota);
    end
        
    OLWH(1,:) =  LWH(1,:) -  margin(1,: )  - margin(2,: ); %��ȣ����ң�
    OLWH(2,:) =  LWH(2,:) -  margin(3,: ) - margin(4,: ); %���ȣ����£�
        
end


% function LWH = LWHunbuffer(LWH, margin)
%     
%     narginchk(2,2);
%     
%     LWH(1,:) =  LWH(1,:) -  margin(1,: )  - margin(2,: ); %��ȣ����ң�
%     LWH(2,:) =  LWH(2,:) -  margin(3,: ) - margin(4,: ); %���ȣ����£�
%     
% end