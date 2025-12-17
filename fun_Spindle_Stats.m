function Output = fun_Spindle_Stats(EEGValues,EEGTimes,EEGArtifact_yn, Fs, MinLength, LT, UT)
% Output = Spindle_Stats(EEGValues,EEGTimes,EEGArtifact_yn, Fs, MinLength, LT, UT)
%
% Script that takes an EEG signal and calculates several sleep spindle
% related numbers.
%
% Dependencies:
% SpindleFilter.m            Filters EEG to spindle frequency
% Logical_Consecutive.m      Determines the length of sleep spindles
% Matlab Signal Processing Toolbox
%
% Inputs:
% EEGValues     1xn Vector      The value of the EEG waveform
% EEGTimes      1xn vector      The time of each EEG datapoint
% Fs            Number          The sampling rate of the EEG waveform
% MinLength     Number          The minimum length for a sleep spindle
%                                  in seconds
% LT            Number          The lower amplitude threshold for detection
%                               (ex 2x median amplitude LT=2)
% UT            Number          The upper amplitude threshold for detection
%                               (ex: toddlers McClain in 6x, adults 9x)
% 
% Outputs:
% Output          A structure with the following fields:
% .Values         The filtered EEG waveform
% .Times          The Times of the filtered EEG waveform
% .Amplitude      The amplitude of the resulting waveform
% .Phase          The phase of the resulting waveform
% .Threshold      The threshold used on the Amplitude to determine a 
%                    sleep spindle.
% .SpindleLogical A logical that reports whether each data point is during
%                    a sleep spindle
% .SpindleIndex   The same as SpindleLogical, except instead of all
%                    spindles being 1, they are an integer of what number
%                    sleep spindle it is. Helps for indexing particular
%                    sleep spindles.
% 
% Given K sleep spindles detected
% 
% .SB.StartTime   1xK vector   The start time of the sleep spindle.
% .SB.EndTime     1xK vector   The end time of the sleep spindle.
% .SB.Duration    1xK vector   The duration of the sleep spindle.
% .SB.Amp         1xK cell     The amplitude of the sleep spindle at each
%                                 datapoint.
% .SB.Phase       1xK cell     The phase of the sleep spindle at each
%                                 datapoint.
% .SB.medAmp      1xK vector   The median amplitude of the sleep spindle.
% .SB.maxAmp      1xK vector   The maximum amplitude of the sleep spindle.
% .SB.sumPhase    1xK vector   The sum of all phase changes for the spindle
%                                 burst.
% .SB.nCycles     1xK vector   The number of cycles for each sleep spindle.
%
% Initally used in: 
% % Sokoloff G, Dooley JC, Glanz RM, Yen RY, Hickerson MM, Evans L, Laughlin
%    HM, Apfelbaum KS, and Blumberg MS (2021). Twitches emerge postnatally 
%    during quiet sleep in human infants and are synchronized with sleep 
%    spindles. Current Biology. https://doi.org/10.1016/j.cub.2021.05.038
%
% Last updated 6/16/2024 by Taylor Christiansen added custom ampltiude and
% spindle length inputs and expanded frequency from 12-14 to 11-15 Hz
%
% Used in Christiansen et al., 2026



[fValues,fTimes] = SpindleFilter(EEGValues,EEGTimes,Fs); % Filter the signal
fValuesAmp = abs(hilbert(fValues)); % Amplitude
fValuesPhase = angle(hilbert(fValues)); % Phase
dPhase = diff(unwrap(fValuesPhase)); % difference in phase
dPhase = [0 dPhase]; % Pad the vector

%% Only include signal not during noise to calculate the median Amplitude
for t = 1:length(fTimes) %loop through EEG data row
    Amp = fValuesAmp (t); % Amplitude at time t
    Art = EEGArtifact_yn (t); % Art y or n at time t
        if Art == 1 ; %if the amplitude is during artifact period
            fValuesAmp_clean(t) = NaN;
        else
            fValuesAmp_clean(t) = Amp;
        end
end

medAmp = nanmedian(fValuesAmp_clean);
LowerThreshold = medAmp*LT; % Change to change threshold
UpperThreshold = medAmp*UT;

fValuesAmpLogical = fValuesAmp_clean > LowerThreshold;
ThreshDuration = Logical_Consecutive(fValuesAmpLogical);
MinDuration = round(Fs*MinLength);

DurationLogical = ThreshDuration.ConsecutiveOutput < MinDuration;
fValuesAmpLogical(DurationLogical & fValuesAmpLogical) = 0;

ThreshDuration = Logical_Consecutive(fValuesAmpLogical);
DurationLogical = ThreshDuration.ConsecutiveOutput < MinDuration;
fValuesAmpLogical(DurationLogical & ~fValuesAmpLogical) = 1;

ThreshDuration = Logical_Consecutive(fValuesAmpLogical);

dfValuesAmpLogical = diff(fValuesAmpLogical);
SSI_raw = find(dfValuesAmpLogical == 1); % Spindle Start Index
SEI_raw = find(dfValuesAmpLogical == -1); % Spindle End Index

nSpindles_raw = min([length(SEI_raw); length(SSI_raw)]);

SSI = []; %initalize
SEI = [];

%Loop Through existing spindles
for iSpindle = 1:nSpindles_raw
Amp_loop = [];
Amp_loop = fValuesAmp(SSI_raw(iSpindle):SEI_raw(iSpindle));
st = SSI_raw (iSpindle);
en = SEI_raw (iSpindle);
% Check if any value in the event crosses the Upper Threshold
    if any(Amp_loop > UpperThreshold) % if at any point the amplitude of spindle crosses the UpperThreshold save it
        % If it does, include the entire event
        SSI = [SSI st];
        SEI = [SEI en];
    end
end

nSpindles = [];
nSpindles = min([length(SEI); length(SSI)]);

Output.Values = fValues;
Output.Times = fTimes;
Output.Amplitude = fValuesAmp_clean;
Output.Phase = fValuesPhase;
Output.LowerThreshold = LowerThreshold;
Output.UpperThreshold = UpperThreshold;
Output.SpindleLogical = fValuesAmpLogical;
Output.SpindleIndex = zeros(size(fValuesAmpLogical));

if nSpindles > 0
Output.nSpindles = nSpindles;
else
     Output.nSpindles = 0;
end

% Save the information for the identified sleep spindles
for iSpindle = 1:nSpindles
    Output.SpindleIndex(SSI(iSpindle):SEI(iSpindle)) = iSpindle; % Create string for easy sleep spindle indexing
    Output.SB.StartTime(iSpindle) = SSI(iSpindle)/Fs; % Get sleep spindle start time (in seconds)
    Output.SB.EndTime(iSpindle) = SEI(iSpindle)/Fs; % Get sleep spindle end time (in seconds)
    Output.SB.Duration(iSpindle) = Output.SB.EndTime(iSpindle) - Output.SB.StartTime(iSpindle); % Get sleep spindle duration
    Output.SB.Amp{iSpindle} = fValuesAmp(SSI(iSpindle):SEI(iSpindle)); % Get sleep spindle amplitude throughout 
    Output.SB.Phase{iSpindle} = fValuesPhase(SSI(iSpindle):SEI(iSpindle)); % Get sleep spindle phase throughout
    Output.SB.medAmp(iSpindle) = median(fValuesAmp(SSI(iSpindle):SEI(iSpindle))); % Get sleep spindle median amplitude
    Output.SB.maxAmp(iSpindle) = max(fValuesAmp(SSI(iSpindle):SEI(iSpindle))); % Get sleep spindle max amplitude
    Output.SB.sumPhase(iSpindle) = sum(dPhase(SSI(iSpindle):SEI(iSpindle))); % Get the phase duration of the sleep spindle
    Output.SB.nCycles(iSpindle) = Output.SB.sumPhase(iSpindle) / (2*pi); % Get the number of cycles of the sleep spindle
    Output.SB.frequencies(iSpindle)= Output.SB.nCycles(iSpindle)/Output.SB.Duration(iSpindle); %use the number of cycles and the duration of the spindle to calculate frequency
end

end




