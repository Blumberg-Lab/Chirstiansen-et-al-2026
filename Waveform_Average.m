function WaveformEvent = Waveform_Average(EventTimes,WaveTimes,WaveValues,Window,PreWindow)
% function WaveformEvent = Waveform_Average(EventTimes,WaveTimes,WaveValues,Window,PreWindow)
%
% This is a function that extracts event-triggered waveform segments from a
% large waveform vector, and outputs times (relative to the event), trial
% by trial values, and mean overall value.
% 
% Input: EventTimes         Time of triggering stimulus (seconds)
%        WaveTimes          Times of waveform datapoints (seconds)
%        WaveValues         Values of waveform datapoints (seconds)
%        Window             Data window for waveform averages (seconds)
%        PreWindow          Time before trigger in window (seconds)
%
% Output: A structure with the following values
%    .MTX                   A matrix (nEventTimes by nWindowBins) of trial
%                           by trial waveform values.
%    .WaveformAverage       The mean waveform across all trials
%    .WindowBins            The bins of each datapoint (seconds)
%
% This script can be used (with caution) for waveforms with uneven sampling
% rates, as it determines the average period between datapoints
% when calculating the window.
%
% If waveform data isn't available due to event times, window size, and
% available waveform data, it replaces unavailable data with NaN's,
% although this has not yet been tested and may have some off-by-one
% errors.
%
% Written by James Dooley on 6/16/2019.
% Version 1.0


if length(WaveTimes) > 10000
    WavePeriod = mean(diff(WaveTimes(1:10000)));
else
    WavePeriod = mean(diff(WaveTimes));
end

WindowBins = -PreWindow:WavePeriod:(Window-PreWindow);
PreWindowBins = length(WavePeriod:WavePeriod:PreWindow);
PostWindowBins = length(0:WavePeriod:(Window-PreWindow));

WaveformEvent.MTX = zeros(length(EventTimes),length(WindowBins));

for iEvent = 1:length(EventTimes)
    cEvent = EventTimes(iEvent);
    [~, cEventIndex] = min(abs(WaveTimes-cEvent));
    cStartIndex = cEventIndex - PreWindowBins;
    cEndIndex = cStartIndex+(length(WindowBins)-1);
    
    if cStartIndex > 0 && cEndIndex < length(WaveValues)
        WaveformEvent.MTX(iEvent,:) = WaveValues(cStartIndex:cEndIndex);
    elseif cStartIndex < 0
        nNaNs = abs(cStartIndex)+1;
        WaveformEvent.MTX(iEvent,1:nNaNs) = NaN;
        WaveformEvent.MTX(iEvent,nNaNs+1:end) = WaveValues(1:cEndIndex);
    elseif cEndIndex > length(WaveValues)
        nNaNs = cEndIndex - length(WaveValues);
        WaveformEvent.MTX(iEvent,(end-(nNaNs+1)):end) = NaN;
        WaveformEvent.MTX(iEvent,1:end-nNaNs) = WaveValues(cStartIndex:length(WaveValues));
    end
end

WaveformEvent.WindowBins = WindowBins;
WaveformEvent.WaveformAverage = nanmean(WaveformEvent.MTX);