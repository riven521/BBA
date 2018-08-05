%%PARAMETERNITIALIZE Initialize the parameter data structure.
%
%% Form
%  p = ParameterInitialize( varargin )
%
%% Description
% Initializes the algorithm's parameter data structure using parameter pairs.
%
%% Inputs
% varargin:  ('parameter',value,...)
%
% 'whichStripH'                             (1,1) 
% 'whichBinH'                               (1,1) 
% 'whichSortItemOrder'               (1,1)
% 'whichRotation'                         (1,1)
% 'whichRotationHori'                  (1,1)
%
%% Outputs
%   p	(.)  Data structure
function p = ParameterInitialize( varargin )

% Defaults
p.whichStripH = 1;
p.whichBinH = 1;
p.whichSortItemOrder = 1;
p.whichRotation = 1;
p.whichRotationHori = 1;

n = length(varargin);

for k = 1:2:length(varargin)
    switch varargin{k}
        case 'whichStripH'
            p.whichStripH                        = varargin{k+1};
        case 'whichBinH'
            p.whichBinH                         = varargin{k+1};
        case 'whichSortItemOrders'
            p.whichSortItemOrder          = varargin{k+1};
        case 'whichRotation'
            p.whichRotation                 = varargin{k+1};
        case 'whichRotationHori'
            p.whichRotationHori         = varargin{k+1};
    end
end

end