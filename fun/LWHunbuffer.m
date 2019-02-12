function OLWH = LWHunbuffer(LWH, margin, rota)
% remove LWH的buffer 且恢复旋转状态
% 若没有参数3，则不恢复旋转

    narginchk(2,4);
    
    OLWH = LWH;     
    
    if nargin > 2
        rota = logical(rota); 
        LWH(1, rota ) = OLWH(2, rota);
        LWH(2, rota ) = OLWH(1, rota);
    end
        
    OLWH(1,:) =  LWH(1,:) -  margin(1,: )  - margin(2,: ); %宽度（左右）
    OLWH(2,:) =  LWH(2,:) -  margin(3,: ) - margin(4,: ); %长度（上下）
        
end


% function LWH = LWHunbuffer(LWH, margin)
%     
%     narginchk(2,2);
%     
%     LWH(1,:) =  LWH(1,:) -  margin(1,: )  - margin(2,: ); %宽度（左右）
%     LWH(2,:) =  LWH(2,:) -  margin(3,: ) - margin(4,: ); %长度（上下）
%     
% end