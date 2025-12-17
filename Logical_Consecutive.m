function Data = Logical_Consecutive(LogicalInput, varargin)
% Data = Logical_Consecutive(LogicalInput, varargin)
%
% This takes a logical input (only 1s and 0s) and gives you the value of
% the first data point (1 or 0), the number of repeats at each value, and a
% modified vector that replaces 1s and 0s with the duration of that value.
% 
% Input:    LogicalInput      A 1xn vector of 1s and 0s
% Optional: EndsNaN           0 by default, replaces the ending values with
%                                 NaN
% 
% Output: Data                A structure with the following fields:
%
%         .FirstDataPoint     The value of the first number (0 or 1)
%         .ConsecutiveOutput  A 1xn string with each value replaced with
%                                an integer of the number of repeats
%         .nRepeats           A vector of the length of each repeat.
%
%
% Example: 
% LogicalInput =  [0 0 0 1 1 1 1 0 0 0 0 0 1 1 1 0 0 1]
%
% Data = Logical_Consecutive(LogicalInput)
% Data.FirstDataPoint =    [0]
% Data.ConsecutiveOutput = [3 3 3 4 4 4 4 5 5 5 5 5 3 3 3 2 2 1]
% Data.nRepeats = [3 4 5 3 2 1]
%
% Data = Logical_Consecutive(LogicalInput,'EndsNaN',1)
% Data.FirstDataPoint =    [0]
% Data.ConsecutiveOutput = [NaN NaN NaN 4 4 4 4 5 5 5 5 5 3 3 3 2 2 NaN]
% Data.nRepeats = [NaN 4 5 3 2 NaN]
%
% Modified from Dooley, Glanz, et al. 2020
% Last updated 2/23/2021 by Jimmy Dooley (james-c-dooley@uiowa.edu)
%
% CC-by-4.0 to Jimmy Dooley
% Used in Sokoloff et al., 2021.
%
% Cite as:
% Sokoloff G, Dooley JC, Glanz RM, Yen RY, Hickerson MM, Evans L, Laughlin
%    HM, Apfelbaum KS, and Blumberg MS (2021). Twitches emerge postnatally 
%    during quiet sleep in human infants and are synchronized with sleep 
%    spindles. Current Biology. https://doi.org/10.1016/j.cub.2021.05.038

Params = inputParser;
Params.addRequired('LogicalInput', @(x) islogical(x) | isnumeric(x));
Params.addParameter('EndsNaN', 0, @isnumeric);
Params.parse(LogicalInput, varargin{:});

EndsNaN = Params.Results.EndsNaN;

LogicalInput = logical(LogicalInput);

FirstDataPoint = LogicalInput(1);
dLogicalInput = abs(diff(LogicalInput));
dLogicalInputIndex = find(dLogicalInput == 1);
if sum(dLogicalInput) > 0
    dLogicalInputIndex(end+1) = (length(LogicalInput) - max(dLogicalInput))+1;
end
if ~isempty(dLogicalInputIndex)
    d_dLogicalInputIndex = diff(dLogicalInputIndex);

    Output = zeros(size(LogicalInput));
    if ~isempty(d_dLogicalInputIndex)
        Output(1:dLogicalInputIndex) = dLogicalInputIndex(1);
        for iIndex = 1:length(dLogicalInputIndex)-1
            Output((dLogicalInputIndex(iIndex)+1):dLogicalInputIndex(iIndex+1)) = d_dLogicalInputIndex(iIndex);
        end
        Output(dLogicalInputIndex(end):end) = d_dLogicalInputIndex(end);

        Data.FirstDataPoint = FirstDataPoint;
        Data.ConsecutiveOutput = Output;
        Data.nRepeats = [Data.ConsecutiveOutput(1), d_dLogicalInputIndex];
    else

        Data.FirstDataPoint = FirstDataPoint;
        Data.ConsecutiveOutput = ones(size(LogicalInput));
        Data.ConsecutiveOutput(1:dLogicalInputIndex) = dLogicalInputIndex;
        Data.ConsecutiveOutput(dLogicalInputIndex+1:end) = length(LogicalInput)-dLogicalInputIndex;        
        Data.nRepeats = [Data.ConsecutiveOutput(1), Data.ConsecutiveOutput(end)];
    end
else
    
    Data.FirstDataPoint = FirstDataPoint;
    Data.ConsecutiveOutput = ones(size(LogicalInput))*length(LogicalInput);
    Data.nRepeats = length(LogicalInput);
end

if EndsNaN
    Data.ConsecutiveOutput(1:Data.ConsecutiveOutput(1)) = NaN;
    if ~isempty(dLogicalInputIndex)
        Data.ConsecutiveOutput(end-Data.ConsecutiveOutput(end):end) = NaN;
    end
    Data.nRepeats([1 end]) = NaN;
end
