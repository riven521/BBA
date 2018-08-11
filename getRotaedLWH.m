% % V2 : margin: ÿ�����̵�margin . ������߶Ե�,����Rotaed����; 
function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, margin)
        % ����ת�Ĳ���������buff������ת�ĳ����滻������buff��
        LWHRotaed = LWH;
        LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - (margin(3, Rotaed) +margin(4, Rotaed)) + (margin(1, Rotaed) +margin(2, Rotaed));
        LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - (margin(1, Rotaed) +margin(2, Rotaed)) + (margin(3, Rotaed) +margin(4, Rotaed)); 
end

% % V1 : buff: ���̼�ļ�϶ . ������߶Ե�,����Rotaed����; 
% function [LWHRotaed]  =getRotaedLWH(LWH,Rotaed, Buff)
%         % ����ת�Ĳ���������buff������ת�ĳ����滻������buff��
%         if nargin < 3
%             Buff = zeros(size(LWH));
%         end
%         LWHRotaed = LWH;
%         LWHRotaed(1, Rotaed) = LWH(2, Rotaed) - Buff(2, Rotaed) +Buff(1, Rotaed);
%         LWHRotaed(2, Rotaed) = LWH(1, Rotaed) - Buff(1, Rotaed) +Buff(2, Rotaed);      
% end