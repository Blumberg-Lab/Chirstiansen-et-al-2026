function [pSpindle] = fun_pSpindle(time,twitchforspindle, SpindleFinal, chan)
% fun_pSpindle finds the probability that a spindles will occur in relation
% to the onset of a twitch
%
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
%   chan= chan number that you want to calculate pspindle for (ex: 36 is
%   channel 36)
%   
%   Output: 
%   is a structure PSpindle 
%   pSpindle.(NREM or N2 or N3).Twitchwithspindle.withspindle_number =  number of
%   twitches which co-occur with a spindle
%   pSpindle.(NREM or N2 or N3).Twitchwithspindle.totalnumber = total number of
%   twitches in that state
%   pSpindle.(NREM or N2 or N3).Twitchwithspindle.withspindle_proportion =
%   PSpindle|Twitch
%   pSpindle.(NREM or N2 or N3).Twitchwithspindle.withspindle_on = list of twitches
%   onsets for twitches that co-occur with spindles
%   pSpindle.(NREM or N2 or N3).Twitchwithspindle.nospindle_on = list of twitches
%   onsets for twitches that do not co-occur with spindles
%   pSpindle.(NREM or N2 or N3).WaveformAverage = the pSpindle in +/- 30 sec 
%   from the onset of a twitch across time
%
%   Used in Christiansen et al., 2026


sessionduration = time(end);

cLength = round(sessionduration*10); 
Times = 0.1:0.1:cLength/10;
curChan=['E' num2str(chan)];
thisSpindle=SpindleFinal.(curChan).bystate;

for st = 1:3

if st ==1 

sp_state = ['N' num2str(2)]; 

cSpindles = thisSpindle.(sp_state);
nSpindles=length(cSpindles.StartTime); % Get the number of spindles
SpindleLogical = zeros(1,cLength);

for iSpindle = 1:nSpindles % Loop through all the spindles
    cStart = round(cSpindles.StartTime(iSpindle)*10); % Index for this spindle's start time
    cEnd = round(cSpindles.EndTime(iSpindle)*10); % Index for this spindle's end time
    SpindleLogical(cStart:cEnd) = 1; % For that spindle, set the values when it's happening to be 1
end

cTwitches = twitchforspindle.bystate.(sp_state); %if you want to do probability for a specific state
nTwitches = length(cTwitches.onset); % Get the number of spindles
tdur=0.5; %0.5sec

pSpindle.(sp_state) = Waveform_Average(cTwitches.onset,Times,SpindleLogical,40,20);
pSpindle.(sp_state).WaveformAverage_smooth=smoothdata(pSpindle.(sp_state).WaveformAverage,"movmean", 5);

%Initalize;
twitchlog=[];
twitchnolog=[];

for tw=1:length(cTwitches.onset)
    tw_on=cTwitches.onset(tw);
    sp_log=zeros(1, nSpindles);
    for sp=1:nSpindles
        spin_on=cSpindles.StartTime(sp);
        spin_off=cSpindles.EndTime(sp);
        if spin_on<=tw_on && tw_on<=spin_off
            sp_log(1,sp)=1;
        end
    end

    if any(sp_log) == 1
        twitchlog = [twitchlog tw_on];
    else
        twitchnolog = [twitchnolog tw_on];
    end

end

pSpindle.(sp_state).Twitchwithspindle.withspindle_number = length(twitchlog);
pSpindle.(sp_state).Twitchwithspindle.totalnumber = length (cTwitches.onset);
pSpindle.(sp_state).Twitchwithspindle.withspindle_proportion = length(twitchlog)/length (cTwitches.onset); %AKA PSpindle given tw
pSpindle.(sp_state).Twitchwithspindle.withspindle_on = twitchlog;
pSpindle.(sp_state).Twitchwithspindle.nospindle_on = twitchnolog;       
end

if st ==2
sp_state = ['N' num2str(3)]; 
cSpindles = thisSpindle.(sp_state);
nSpindles=length(cSpindles.StartTime); % Get the number of spindles
SpindleLogical = zeros(1,cLength);

for iSpindle = 1:nSpindles % Loop through all the spindles
    cStart = round(cSpindles.StartTime(iSpindle)*10); % Index for this spindle's start time
    cEnd = round(cSpindles.EndTime(iSpindle)*10); % Index for this spindle's end time
    SpindleLogical(cStart:cEnd) = 1; % For that spindle, set the values when it's happening to be 1
end

cTwitches = twitchforspindle.bystate.(sp_state); %if you want to do probability for a specific state
nTwitches = length(cTwitches.onset); % Get the number of spindles
tdur=0.5; %0.5sec

pSpindle.(sp_state) = Waveform_Average(cTwitches.onset,Times,SpindleLogical,40,20);
pSpindle.(sp_state).WaveformAverage_smooth=smoothdata(pSpindle.(sp_state).WaveformAverage,"movmean", 5);

%Initalize;
twitchlog=[];
twitchnolog=[];

for tw=1:length(cTwitches.onset)
    tw_on=cTwitches.onset(tw);
    sp_log=zeros(1, nSpindles);
    for sp=1:nSpindles
        spin_on=cSpindles.StartTime(sp);
        spin_off=cSpindles.EndTime(sp);
        if spin_on<=tw_on && tw_on<=spin_off
            sp_log(1,sp)=1;
        end
    end

    if any(sp_log) == 1
        twitchlog = [twitchlog tw_on];
    else
        twitchnolog = [twitchnolog tw_on];
    end

end

pSpindle.(sp_state).Twitchwithspindle.withspindle_number = length(twitchlog);
pSpindle.(sp_state).Twitchwithspindle.totalnumber = length (cTwitches.onset);
pSpindle.(sp_state).Twitchwithspindle.withspindle_proportion = length(twitchlog)/length (cTwitches.onset); %AKA PSpindle given tw
pSpindle.(sp_state).Twitchwithspindle.withspindle_on = twitchlog;
pSpindle.(sp_state).Twitchwithspindle.nospindle_on = twitchnolog; 
end

if st==3
sp_state = ['NREM']; 
cSpindles = thisSpindle.(sp_state);
nSpindles=length(cSpindles.StartTime); % Get the number of spindles
SpindleLogical = zeros(1,cLength);

for iSpindle = 1:nSpindles % Loop through all the spindles
    cStart = round(cSpindles.StartTime(iSpindle)*10); % Index for this spindle's start time
    cEnd = round(cSpindles.EndTime(iSpindle)*10); % Index for this spindle's end time
    SpindleLogical(cStart:cEnd) = 1; % For that spindle, set the values when it's happening to be 1
end

cTwitches = twitchforspindle.bystate.(sp_state); %if you want to do probability for a specific state
nTwitches = length(cTwitches.onset); % Get the number of spindles
tdur=0.5; %0.5sec

pSpindle.(sp_state) = Waveform_Average(cTwitches.onset,Times,SpindleLogical,40,20);
pSpindle.(sp_state).WaveformAverage_smooth=smoothdata(pSpindle.(sp_state).WaveformAverage,"movmean", 5);

%Initalize;
twitchlog=[];
twitchnolog=[];

for tw=1:length(cTwitches.onset)
    tw_on=cTwitches.onset(tw);
    sp_log=zeros(1, nSpindles);
    for sp=1:nSpindles
        spin_on=cSpindles.StartTime(sp);
        spin_off=cSpindles.EndTime(sp);
        if spin_on<=tw_on && tw_on<=spin_off
            sp_log(1,sp)=1;
        end
    end

    if any(sp_log) == 1
        twitchlog = [twitchlog tw_on];
    else
        twitchnolog = [twitchnolog tw_on];
    end

end

pSpindle.(sp_state).Twitchwithspindle.withspindle_number = length(twitchlog);
pSpindle.(sp_state).Twitchwithspindle.totalnumber = length (cTwitches.onset);
pSpindle.(sp_state).Twitchwithspindle.withspindle_proportion = length(twitchlog)/length (cTwitches.onset); %AKA PSpindle given tw
pSpindle.(sp_state).Twitchwithspindle.withspindle_on = twitchlog;
pSpindle.(sp_state).Twitchwithspindle.nospindle_on = twitchnolog;  




end

end

end