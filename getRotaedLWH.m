% V2 : buff: ÿ�����̵�margin . ������߶Ե�,����Rotaed����; 
% function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, buff)
%         % ����ת�Ĳ���������buff������ת�ĳ����滻������buff��
%         LWHRotaed = LWH;
%         LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - (buff(3, Rotaed) +buff(4, Rotaed)) + (buff(1, Rotaed) +buff(2, Rotaed));
%         LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - (buff(1, Rotaed) +buff(2, Rotaed)) + (buff(3, Rotaed) +buff(4, Rotaed)); 
% end

% V1 : buff: ���̼�ļ�϶ . ������߶Ե�,����Rotaed����; 
function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, Buff)
        % ����ת�Ĳ���������buff������ת�ĳ����滻������buff��
        if nargin < 3
            Buff = zeros(size(LWH));
        end
        LWHRotaed = LWH;
        LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - Buff(2, Rotaed) +Buff(1, Rotaed);
        LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - Buff(1, Rotaed) +Buff(2, Rotaed); 
      
end