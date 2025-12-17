function [pTwitch] = fun_pTwitch(time,twitchforspindle, SpindleFinal, chan)
% fun_pTwitch finds the probability that a twitch will occur in relation
% to a sleep spindle
%   Dependency: Waveform_Average script
%
%   Inputs:
%   time=vector of times in seconds for each EEG sample in the sleep session
%   twitchforspindle= structure with onset times of twitches during
%                     artifact free periods where sleep spindle were identified
%   SpindleFinal = structure with the onset times of sleep spindles in each
%                   indetified electrode channel
%   within both twitchforspindle and spindlefinal twitches/spindle are
%   separated by state so in this script you can call onset of for instance
%   on N2 twitches (ex twitchforspindle.bystate.N2.onset)
%   chan= chan number that you want to calculate pTwitch for (ex: 36 is
%   channel 36)
%   
%   Output: is a structure PTwitch
%   pTwitch.(NREM or N2 or N3).Spindlewithtwitch.withtwitch_number =  number of
%   twitches which co-occur with a spindle
%   pTpindle.(NREM or N2 or N3).Spindlewithtwitch.totalnumber = total number of
%   spindle in that state
%   pTpindle.(NREM or N2 or N3).Twitchwithspindle.withspindle_proportion =
%   PTwitch|Spindle
%   pTwitch.(NREM or N2 or N3).(SpinMid or SpinOn).WaveformAverage = the
%   pTwitch in +/- 30 sec from the onset or midpoint of a spindle across
%   time
%
%   Used in Christiansen et al., 2026

sessionduration = time(end);

cLength = round(sessionduration*10); 
Times = 0.1:0.1:cLength/10;
curChan=['E' num2str(chan)];
thisSpindle=SpindleFinal.E36.bystate;

% Next, we create a big vector in 1/10th of a second incraments that will
% eventually be a 0 if there is no spindle, and a 1 if there is a spindle.
for st = 1:3

if st ==1 

sp_state = ['N' num2str(2)]; 
TwitchLogical = zeros(1,cLength);
cTwitches = twitchforspindle.bystate.(sp_state); %if you want to do probability for a specific state
nTwitches = length(cTwitches.onset); % Get the number of spindles
tdur=0.5; %0.5sec

for iTwitch = 1:nTwitches % Loop through all the spindles
    cStart = round(cTwitches.onset(iTwitch)*10); % Index for this twitch start time
    cEnd = cStart+ (tdur*10); % Index for this twitch end time (start plus 0.5 secs)
    TwitchLogical(cStart:cEnd) = 1; % For that twitch, set the values when it's happening to be 1
end

% And here's where the magic happens. This creates a "waveform average" -
% which functions to give you the probability because we know for every
% tenth of a second, there either is or isn't a twitch. So when we treat
% the logical twitch vector as a waveform, we end up with the probability.
%
% cSpindles is the time of the Spindles (in seconds)
% Times is the time vector (in 1/10th of a second intervals)
% Twitch Logical is what we made above
% 60 is the size of the window (in seconds)
% 30 is the time before the trigger (twitches) to plot

cSpindles = thisSpindle.(sp_state);
durSpindle = cSpindles.Duration;
midSpindles =cSpindles.StartTime + (durSpindle/2); %midpoint
nSpindles=length(cSpindles.StartTime); % Get the number of spindles

pTwitch.(sp_state).SpinOn = Waveform_Average(cSpindles.StartTime,Times,TwitchLogical,60,30);
pTwitch.(sp_state).SpinMid = Waveform_Average(midSpindles,Times,TwitchLogical,60,30);
pTwitch.(sp_state).SpinOn.WaveformAverage_smooth=smoothdata(pTwitch.(sp_state).SpinOn.WaveformAverage,"movmean", 5)
pTwitch.(sp_state).SpinMid.WaveformAverage_smooth=smoothdata(pTwitch.(sp_state).SpinMid.WaveformAverage,"movmean", 5)

figure
plot(pTwitch.(sp_state).SpinOn.WindowBins,pTwitch.(sp_state).SpinOn.WaveformAverage_smooth)
xlabel('Time in relation to spindle onset(sec)')
ylabel('PTwitch')
title(sp_state)
ylim([0 1])

figure
plot(pTwitch.(sp_state).SpinMid.WindowBins,pTwitch.(sp_state).SpinMid.WaveformAverage_smooth)
xlabel('Time in relation to spindle midpoint(sec)')
ylabel('PTwitch')
title(sp_state)
ylim([0 1])

spinlog=[];
spinnolog = [];
for sp=1: length(cSpindles.StartTime)
    spin_on=cSpindles.StartTime(sp);
    spin_off=cSpindles.EndTime(sp);
    tw_log=zeros(1,length(cTwitches.onset));
    for tw=1:length(cTwitches.onset)
        tw_on=cTwitches.onset(tw);
        if spin_on<=tw_on && tw_on<=spin_off
            tw_log(1,tw)=1;
        end
    end

    if any(tw_log) == 1
        spinlog = [spinlog spin_on];
    else
        spinnolog = [spinnolog spin_on];
    end

end

pTwitch.(sp_state).Spindlewithtwitch.withtwitch_number = length(spinlog);
pTwitch.(sp_state).Spindlewithtwitch.totalnumber = nSpindles;
pTwitch.(sp_state).Spindlewithtwitch.withtwitch_proportion = length(spinlog)/nSpindles; %AKA PTwitch given spin
pTwitch.(sp_state).Spindlewithtwitch.withtwitch_on = spinlog;
pTwitch.(sp_state).Spindlewithtwitch.notwitch_on = spinnolog;       
end

if st ==2
sp_state = ['N' num2str(3)]; 
TwitchLogical = zeros(1,cLength);
cTwitches = twitchforspindle.bystate.(sp_state); %if you want to do probability for a specific state
nTwitches = length(cTwitches.onset); % Get the number of spindles
tdur=0.5; %0.5sec

for iTwitch = 1:nTwitches % Loop through all the spindles
    cStart = round(cTwitches.onset(iTwitch)*10); % Index for this twitch start time
    cEnd = cStart+ (tdur*10); % Index for this twitch end time (start plus 0.5 secs)
    TwitchLogical(cStart:cEnd) = 1; % For that twitch, set the values when it's happening to be 1
end

% And here's where the magic happens. This creates a "waveform average" -
% which functions to give you the probability because we know for every
% tenth of a second, there either is or isn't a twitch. So when we treat
% the logical twitch vector as a waveform, we end up with the probability.
%
% cSpindles is the time of the Spindles (in seconds)
% Times is the time vector (in 1/10th of a second intervals)
% Twitch Logical is what we made above
% 60 is the size of the window (in seconds)
% 30 is the time before the trigger (twitches) to plot

cSpindles = thisSpindle.(sp_state);
durSpindle = cSpindles.Duration;
midSpindles =cSpindles.StartTime + (durSpindle/2); %midpoint
nSpindles=length(cSpindles.StartTime); % Get the number of spindles

pTwitch.(sp_state).SpinOn = Waveform_Average(cSpindles.StartTime,Times,TwitchLogical,60,30);
pTwitch.(sp_state).SpinMid = Waveform_Average(midSpindles,Times,TwitchLogical,60,30);
pTwitch.(sp_state).SpinOn.WaveformAverage_smooth=smoothdata(pTwitch.(sp_state).SpinOn.WaveformAverage,"movmean", 5)
pTwitch.(sp_state).SpinMid.WaveformAverage_smooth=smoothdata(pTwitch.(sp_state).SpinMid.WaveformAverage,"movmean", 5)

figure
plot(pTwitch.(sp_state).SpinOn.WindowBins,pTwitch.(sp_state).SpinOn.WaveformAverage_smooth)
xlabel('Time in relation to spindle onset(sec)')
ylabel('PTwitch')
title(sp_state)
ylim([0 1])

figure
plot(pTwitch.(sp_state).SpinMid.WindowBins,pTwitch.(sp_state).SpinMid.WaveformAverage_smooth)
xlabel('Time in relation to spindle midpoint(sec)')
ylabel('PTwitch')
title(sp_state)
ylim([0 1])

spinlog=[];
spinnolog = [];
for sp=1: length(cSpindles.StartTime)
    spin_on=cSpindles.StartTime(sp);
    spin_off=cSpindles.EndTime(sp);
    tw_log=zeros(1,length(cTwitches.onset));
    for tw=1:length(cTwitches.onset)
        tw_on=cTwitches.onset(tw);
        if spin_on<=tw_on && tw_on<=spin_off
            tw_log(1,tw)=1;
        end
    end

    if any(tw_log) == 1
        spinlog = [spinlog spin_on];
    else
        spinnolog = [spinnolog spin_on];
    end

end

pTwitch.(sp_state).Spindlewithtwitch.withtwitch_number = length(spinlog);
pTwitch.(sp_state).Spindlewithtwitch.totalnumber = nSpindles;
pTwitch.(sp_state).Spindlewithtwitch.withtwitch_proportion = length(spinlog)/nSpindles; %AKA PTwitch given spin
pTwitch.(sp_state).Spindlewithtwitch.withtwitch_on = spinlog;
pTwitch.(sp_state).Spindlewithtwitch.notwitch_on = spinnolog;     
end

if st==3
sp_state = ['NREM']; 
TwitchLogical = zeros(1,cLength);
cTwitches = twitchforspindle.bystate.(sp_state); %if you want to do probability for a specific state
nTwitches = length(cTwitches.onset); % Get the number of spindles
tdur=0.5; %0.5sec

for iTwitch = 1:nTwitches % Loop through all the spindles
    cStart = round(cTwitches.onset(iTwitch)*10); % Index for this twitch start time
    cEnd = cStart+ (tdur*10); % Index for this twitch end time (start plus 0.5 secs)
    TwitchLogical(cStart:cEnd) = 1; % For that twitch, set the values when it's happening to be 1
end

% And here's where the magic happens. This creates a "waveform average" -
% which functions to give you the probability because we know for every
% tenth of a second, there either is or isn't a twitch. So when we treat
% the logical twitch vector as a waveform, we end up with the probability.
%
% cSpindles is the time of the Spindles (in seconds)
% Times is the time vector (in 1/10th of a second intervals)
% Twitch Logical is what we made above
% 60 is the size of the window (in seconds)
% 30 is the time before the trigger (twitches) to plot

cSpindles = thisSpindle.(sp_state);
durSpindle = cSpindles.Duration;
midSpindles =cSpindles.StartTime + (durSpindle/2); %midpoint
nSpindles=length(cSpindles.StartTime); % Get the number of spindles

pTwitch.(sp_state).SpinOn = Waveform_Average(cSpindles.StartTime,Times,TwitchLogical,60,30);
pTwitch.(sp_state).SpinMid = Waveform_Average(midSpindles,Times,TwitchLogical,60,30);
pTwitch.(sp_state).SpinOn.WaveformAverage_smooth=smoothdata(pTwitch.(sp_state).SpinOn.WaveformAverage,"movmean", 5)
pTwitch.(sp_state).SpinMid.WaveformAverage_smooth=smoothdata(pTwitch.(sp_state).SpinMid.WaveformAverage,"movmean", 5)

figure
plot(pTwitch.(sp_state).SpinOn.WindowBins,pTwitch.(sp_state).SpinOn.WaveformAverage_smooth)
xlabel('Time in relation to spindle onset(sec)')
ylabel('PTwitch')
title(sp_state)
ylim([0 1])

figure
plot(pTwitch.(sp_state).SpinMid.WindowBins,pTwitch.(sp_state).SpinMid.WaveformAverage_smooth)
xlabel('Time in relation to spindle midpoint(sec)')
ylabel('PTwitch')
title(sp_state)
ylim([0 1])

spinlog=[];
spinnolog = [];
for sp=1: length(cSpindles.StartTime)
    spin_on=cSpindles.StartTime(sp);
    spin_off=cSpindles.EndTime(sp);
    tw_log=zeros(1,length(cTwitches.onset));
    for tw=1:length(cTwitches.onset)
        tw_on=cTwitches.onset(tw);
        if spin_on<=tw_on && tw_on<=spin_off
            tw_log(1,tw)=1;
        end
    end

    if any(tw_log) == 1
        spinlog = [spinlog spin_on];
    else
        spinnolog = [spinnolog spin_on];
    end

end

pTwitch.(sp_state).Spindlewithtwitch.withtwitch_number = length(spinlog);
pTwitch.(sp_state).Spindlewithtwitch.totalnumber = nSpindles;
pTwitch.(sp_state).Spindlewithtwitch.withtwitch_proportion = length(spinlog)/nSpindles; %AKA PTwitch given spin
pTwitch.(sp_state).Spindlewithtwitch.withtwitch_on = spinlog;
pTwitch.(sp_state).Spindlewithtwitch.notwitch_on = spinnolog;     
end

end

end