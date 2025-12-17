function [discrete_state, con_state, sleep_state, twitch, remraw, rem, ITI, Summary] = fun_twitchstateremprocess(FinState, EEG, Twitches, REM, Sleep, VideoShift)
%Processes Raw Twitch and State Data
%   1. Loads in Data from tables for FinState, Twitches, Sleep and matrix
%   VideoShift
%   2. Stores state data in discrete_state structure
%   3. Creates a continous vector of sleep state values to make
%   hypnograms: con_state
%   4. Shifts twitch data based on the lag between video start times and EEG
%   5. Exclude twitches during non sleep or non scoreable periods
%   6. Assigns each twitch to a state and combines all twitches in a state
%       All twitch data is stored in sturcture twitch
%   7. Calculates ITIs (stored in ITI structure)
%   8. Twitch Rate and State Durations are calculated and stored in the
%   structure Summary 
%   
%   Used in Christiansen et al., 2026
%% Load Data
%StateData (FinState) to discrete_state structure
discrete_state.start_t=table2array(FinState(:, 1));
discrete_state.end_t=table2array(FinState(:, 2));
discrete_state.state=table2array(FinState(:, 4));
discrete_state.end_t(length(discrete_state.end_t))= EEG.times_sec(end);
discrete_state.dur = discrete_state.end_t - discrete_state.start_t
num = length(discrete_state.start_t);

%TwitchData (Twitches) to twitchraw sturcture
twitchraw.bodypart = table2array(Twitches(:, 1));
twitchraw.side = table2array(Twitches(:, 2));
twitchraw.uaonset = table2array(Twitches(:, 3)); %unadjusted onset
twitchraw.uaoffset = table2array(Twitches(:, 4)); %unadjusted offset
twitchraw.dur = table2array(Twitches(:, 5));
twitchraw.total = length(twitchraw.uaonset);

%REMData (REM) to remraw structure
remraw.uaonset = table2array(REM(:, 1));%unadjusted onset
remraw.uaoffset = table2array(REM(:, 2));%unadjusted offset
remraw.dur = table2array(REM(:, 3));
remraw.total = length(remraw.uaonset);

%SleepWake Data (Sleep) to sleep structure
sleep_state.start_t = table2array(Sleep(:, 1)); %state onset times
sleep_state.end_t = table2array(Sleep(:, 2)); %state offset times
sleep_state.yn = table2array(Sleep(:, 4)); %sleep yn
sleep_state.end_t(length(sleep_state.end_t))= EEG.times_sec(end);
sleep_num = length(sleep_state.start_t);

%Load Video Shift Time Data
shift.vstart_raw = VideoShift(:, 1); %video start times (in Netstation time)
shift.vend_raw = VideoShift(:, 2); %video end times (in Netstation time)
shift.vstart = (shift.vstart_raw-1000)/1000000; %convert netstation time to sec
shift.vend = [];
    
%Replace NaN in Video End with Blank 
for i= 1:length(shift.vend_raw)
    v=shift.vend_raw(i)
if isnan(sym(v))== 1
else
    shift.vend(i) = (v+1000)/1000000; %convert netstation time to sec
end
end
shift.vend = shift.vend'

%% Create Continous State from discrete state (for hypnograms)

%Create Variables
e = discrete_state.end_t(end, :); %end time of recording
con_state.time = 0:0.01:e; %time from start to end in steps of 0.01 sec
state_l =length(con_state.time);
num = length(discrete_state.start_t);

%Initalize variables
con_state.hypno =zeros(1, state_l);  

%Loops through continous time when time falls inbetween a state range the
%variable con_state.hypno assigned the value of that state
for time_row= 1:state_l
t=con_state.time(:, time_row); 
    for on_row= 1:num
     on = discrete_state.start_t(on_row, :); %looping through each state change
     off = discrete_state.end_t(on_row, :);
     sleep = discrete_state.state(on_row, :);
     if on<=t %time is between the onset and offset time of a state
         if t<=off
             con_state.hypno(:,time_row)= sleep; %assign hypno the categorical value of the state 
         else
         end %2nd if 
     else 
     end %1st if
     end %on row
 end %time row
 
 %% Shift REM times
 %Adds video start lag to REMs times to align twitches with EEG
for r=1:remraw.total %loops through twitches
    vsta=shift.vstart(1); %lag time of start of first video
    re_on=remraw.uaonset(r); %original twitch onset
    re_off=remraw.uaoffset(r);%original twitch offset
    re_on_adj = (re_on + vsta);% adds video lag time to twitch onset
    re_off_adj = (re_off + vsta);
    if length(shift.vstart)>1 % if there was more than one video start
        for i=1:length(shift.vstart)-1 %loop through each video -1
            vs = shift.vstart(i+1);
            ve = shift.vend(i);
            pause_curr = vs-ve; %finds the difference between end of video and start of next video
            if re_on > vs % if a twitch is after the new video start the additional lags will be added
            re_on_adj = (re_on_adj + pause_curr); %continues to add video lag time for each additional video in loop i
            re_off_adj = (re_off_adj + pause_curr);
            else
            end 
        end
    else

    end
    remraw.onset(r, 1)=re_on_adj; %saves the adjusted onset
    remraw.offset(r, 1)=re_off_adj;
end

%% Shift Twitch Times
%Adds video start lag to twitch times to align twitches with EEG
for t=1:twitchraw.total %loops through twitches
    vsta=shift.vstart(1); %lag time of start of first video
    tw_on=twitchraw.uaonset(t); %original twitch onset
    tw_off=twitchraw.uaoffset(t);%original twitch offset
    tw_on_adj = (tw_on + vsta);% adds video lag time to twitch onset
    tw_off_adj = (tw_off + vsta);
    if length(shift.vstart)>1 % if there was more than one video start
        for i=1:length(shift.vstart)-1 %loop through each video -1
            vs = shift.vstart(i+1);
            ve = shift.vend(i);
            pause_curr = vs-ve; %finds the difference between end of video and start of next video
            if tw_on > vs % if a twitch is after the new video start the additional lags will be added
            tw_on_adj = (tw_on_adj + pause_curr); %continues to add video lag time for each additional video in loop i
            tw_off_adj = (tw_off_adj + pause_curr);
            else
            end 
        end
    else

    end
    twitchraw.onset(t, 1)=tw_on_adj; %saves the adjusted onset
    twitchraw.offset(t, 1)=tw_off_adj;
end
 
%% Exclude REMs during periods of wake or artifact
%Initalize Variables
remraw.state=[];
remraw.ordinal=[];
remraw.sleepyn_state = [];
remraw.sleepyn_ordinal = [];

rem =[];

rem.onset = [];
rem.offset = [];
rem.dur = [];
rem.state = [];
rem.ordinal = [];

for r = 1:remraw.total %looping through each rem
ron=remraw.onset(r); 
roff=remraw.offset(r);
rdur=remraw.dur(r);

    for row = 1:sleep_num
     on = sleep_state.start_t(row, :); 
     off = sleep_state.end_t(row, :);
     sleepyn = sleep_state.yn(row, :);
     if on<=ron && ron<off %if twitch onset is between the onset and offset time of a state
        remraw.sleepyn_state(r) = sleepyn; %save the sleep_yn for the twitch
        remraw.sleepyn_ordinal(r) = row; %save the sleep_yn row for the twitch
     else 
     end %1st if
     end %row

     for state_row = 1:num
     state_on = discrete_state.start_t(state_row, :); 
     state_off = discrete_state.end_t(state_row, :);
     st = discrete_state.state(state_row, :);
     if state_on<=ron && ron<state_off %if twitch onset is between the onset and offset time of a state
        remraw.state(r) = st; %save the sleep_yn for the twitch
        remraw.ordinal(r) = state_row; %save the sleep_yn row for the twitch
     else 
     end %1st if
     end %row
     
     s=remraw.state(r); % s is the state value of the rem in row r
     o=remraw.ordinal(r); % o is the ordinal state row value of the rem in row r
     sly=remraw.sleepyn_state(r); % sly is the sleep_yn value of the rem in row r
     
     if s == 1 | s == 2 %if rem is during any sleep 
     else
                 % Make new structure twitch that includes only rem during
        % scoreable sleep (structure is rem)
        rem.state=[rem.state s];
        rem.ordinal= [rem.ordinal o];
        rem.onset= [rem.onset ron];
        rem.offset = [rem.offset roff];
        rem.dur = [rem.dur rdur];
     end
 end %t
        rem.state= rem.state'
        rem.ordinal= rem.ordinal';
        rem.onset=rem.onset' ;
        rem.offset = rem.offset';
        rem.dur = rem.dur';
%% Loop through each twitch to assign it a state then use the state data  
%Initalize Variables
twitchraw.state=[];
twitchraw.ordinal=[];
twitchraw.sleepyn_state = [];
twitchraw.sleepyn_ordinal = [];

twitch =[];
twitch.sleepyn_ordinal = [];

twitch.onset = [];
twitch.offset = [];
twitch.dur = [];
twitch.bodypart = [];
twitch.side = [];
twitch.state = [];
twitch.ordinal = [];

twitch.bystate=[];

twitch.bystate.REM=[];
twitch.bystate.REM.onset = [];
twitch.bystate.REM.offset = [];
twitch.bystate.REM.dur = [];
twitch.bystate.REM.bodypart = [];
twitch.bystate.REM.side = [];

twitch.bystate.N1=[]
twitch.bystate.N1.onset = [];
twitch.bystate.N1.offset = [];
twitch.bystate.N1.dur = [];
twitch.bystate.N1.bodypart = [];
twitch.bystate.N1.side = [];

twitch.bystate.N2=[]
twitch.bystate.N2.onset = [];
twitch.bystate.N2.offset = [];
twitch.bystate.N2.dur = [];
twitch.bystate.N2.bodypart = [];
twitch.bystate.N2.side = [];

twitch.bystate.N3=[]
twitch.bystate.N3.onset = [];
twitch.bystate.N3.offset = [];
twitch.bystate.N3.dur = [];
twitch.bystate.N3.bodypart = [];
twitch.bystate.N3.side = [];

twitch.bystate.NREM=[];
twitch.bystate.NREM.onset = [];
twitch.bystate.NREM.offset = [];
twitch.bystate.NREM.dur = [];
twitch.bystate.NREM.bodypart = [];
twitch.bystate.NREM.side = [];

%This loop gets rid of any twitches in twitch raw during non-sleep periods
%Then it documents the row of each twitch in the state column (will be
%used later in ITI script). Then it organizes each twitch by state. 

for t = 1:twitchraw.total %looping through each twitch
twon=twitchraw.onset(t); %onset for each loop
twoff=twitchraw.offset(t); %offset for each loop
twdur=twitchraw.dur(t); %duration for each loop
twbp=twitchraw.bodypart(t); %bodypart for each loop
twside=twitchraw.side(t); %side for each loop

%Looping through sleep_yn bs
    for row = 1:sleep_num
     on = sleep_state.start_t(row, :); 
     off = sleep_state.end_t(row, :);
     sleepyn = sleep_state.yn(row, :);
     if on<=twon & twon<off %if twitch onset is between the onset and offset time of a state
        twitchraw.sleepyn_state(t) = sleepyn; %save the sleep_yn for the twitch
        twitchraw.sleepyn_ordinal(t) = row; %save the sleep_yn row for the twitch
     else 
     end %1st if
     end %row
     
%Looping through combined EEG state
    for on_row= 1:num
     on = discrete_state.start_t(on_row, :); %looping through each state change
     off = discrete_state.end_t(on_row, :);
     sleep = discrete_state.state(on_row, :);
     if on<=twon & twon<off %if twitch onset is between the onset and offset time of a state
        twitchraw.state(t)= sleep; %Save state of each twitch
        twitchraw.ordinal(t)= on_row; %Save what state row the twitch is in 
     else 
     end %1st if
     end %on row
     
     s=twitchraw.state(t); % s is the state value of the twitch in row t
     o=twitchraw.ordinal(t); % ordinal for finstate
     
     sly=twitchraw.sleepyn_state(t); % sly is the sleep_yn value of the twitch in row t
     r=twitchraw.sleepyn_ordinal(t); % row for sleep_yn
     
     %Save only twitches during sleep in the final twitch file
     if s == 1 | s == 2 %if twitch is during any sleep 
     else
        % Make new structure twitch that includes only twitches during
        % scoreable sleep (structure is twitch)
        twitch.state=[twitch.state s];
        twitch.ordinal= [twitch.ordinal o];
        twitch.sleepyn_ordinal = [twitch.sleepyn_ordinal r]
        twitch.onset= [twitch.onset twon];
        twitch.offset = [twitch.offset twoff];
        twitch.dur = [twitch.dur twdur];
        twitch.bodypart = [twitch.bodypart twbp];
        twitch.side = [twitch.side twside];
     end
     %Now sort twitches by state and store in twitch.bystate structure
     
     if s == 0 %If twitch is during REM period
         %add twitch to list of all twitches during the state
         twitch.bystate.REM.onset = [twitch.bystate.REM.onset twon];
         twitch.bystate.REM.offset = [twitch.bystate.REM.offset twoff];
         twitch.bystate.REM.dur = [twitch.bystate.REM.dur twdur];
         twitch.bystate.REM.bodypart = [twitch.bystate.REM.bodypart twbp];
         twitch.bystate.REM.side = [twitch.bystate.REM.side twside];
     else
     end
     
     if s == -1 %If twitch is during N1
         twitch.bystate.N1.onset = [twitch.bystate.N1.onset twon];
         twitch.bystate.N1.offset = [twitch.bystate.N1.offset twoff];
         twitch.bystate.N1.dur = [twitch.bystate.N1.dur twdur];
         twitch.bystate.N1.bodypart = [twitch.bystate.N1.bodypart twbp];
         twitch.bystate.N1.side = [twitch.bystate.N1.side twside];
     else
     end
     
     if s == -2 %If twitch is during N2
         twitch.bystate.N2.onset = [twitch.bystate.N2.onset twon];
         twitch.bystate.N2.offset = [twitch.bystate.N2.offset twoff];
         twitch.bystate.N2.dur = [twitch.bystate.N2.dur twdur];
         twitch.bystate.N2.bodypart = [twitch.bystate.N2.bodypart twbp];
         twitch.bystate.N2.side = [twitch.bystate.N2.side twside];
     else
     end
     
     if s == -3 %If twitch is during N3
         twitch.bystate.N3.onset = [twitch.bystate.N3.onset twon];
         twitch.bystate.N3.offset = [twitch.bystate.N3.offset twoff];
         twitch.bystate.N3.dur = [twitch.bystate.N3.dur twdur];
         twitch.bystate.N3.bodypart = [twitch.bystate.N3.bodypart twbp];
         twitch.bystate.N3.side = [twitch.bystate.N3.side twside];
     else
     end
     
     if s == -1 | s == -2 | s == -3 %If twitch is during any NREM
         twitch.bystate.NREM.onset = [twitch.bystate.NREM.onset twon];
         twitch.bystate.NREM.offset = [twitch.bystate.NREM.offset twoff];
         twitch.bystate.NREM.dur = [twitch.bystate.NREM.dur twdur];
         twitch.bystate.NREM.bodypart = [twitch.bystate.NREM.bodypart twbp];
         twitch.bystate.NREM.side = [twitch.bystate.NREM.side twside];
     else
     end
     
 end %t

%Change orientation
twitch.sleepyn_ordinal = twitch.sleepyn_ordinal';

twitch.onset = twitch.onset';
twitch.offset = twitch.offset';
twitch.dur = twitch.dur';
twitch.bodypart = twitch.bodypart';
twitch.side = twitch.side';
twitch.state = twitch.state';
twitch.ordinal = twitch.ordinal';

twitch.bystate.REM.onset = twitch.bystate.REM.onset';
twitch.bystate.REM.offset = twitch.bystate.REM.offset';
twitch.bystate.REM.dur = twitch.bystate.REM.dur';
twitch.bystate.REM.bodypart = twitch.bystate.REM.bodypart';
twitch.bystate.REM.side = twitch.bystate.REM.side';

twitch.bystate.N1.onset = twitch.bystate.N1.onset';
twitch.bystate.N1.offset = twitch.bystate.N1.offset';
twitch.bystate.N1.dur = twitch.bystate.N1.dur';
twitch.bystate.N1.bodypart = twitch.bystate.N1.bodypart';
twitch.bystate.N1.side = twitch.bystate.N1.side';

twitch.bystate.N2.onset = twitch.bystate.N2.onset';
twitch.bystate.N2.offset = twitch.bystate.N2.offset';
twitch.bystate.N2.dur = twitch.bystate.N2.dur';
twitch.bystate.N2.bodypart = twitch.bystate.N2.bodypart';
twitch.bystate.N2.side = twitch.bystate.N2.side';

twitch.bystate.N3.onset = twitch.bystate.N3.onset';
twitch.bystate.N3.offset = twitch.bystate.N3.offset';
twitch.bystate.N3.dur = twitch.bystate.N3.dur';
twitch.bystate.N3.bodypart = twitch.bystate.N3.bodypart';
twitch.bystate.N3.side = twitch.bystate.N3.side';

twitch.bystate.NREM.onset = twitch.bystate.NREM.onset';
twitch.bystate.NREM.offset = twitch.bystate.NREM.offset';
twitch.bystate.NREM.dur = twitch.bystate.NREM.dur';
twitch.bystate.NREM.bodypart = twitch.bystate.NREM.bodypart';
twitch.bystate.NREM.side = twitch.bystate.NREM.side';

 %% Total number of twitches in structure twitch (scored during sleep and
 %visible/good EEG)
 Summary.TwitchTotal.Sleep = length(twitch.onset);
 Summary.TwitchTotal.REM = length(twitch.bystate.REM.onset);
 Summary.TwitchTotal.N1 = length(twitch.bystate.N1.onset);
 Summary.TwitchTotal.N2 = length(twitch.bystate.N2.onset);
 Summary.TwitchTotal.N3 = length(twitch.bystate.N3.onset);
 Summary.TwitchTotal.NREM = length(twitch.bystate.NREM.onset);
 
 %% Seperate twitch by sleepyn ordinal and then calculate ITIs
 %This allows ITIs between sleep wake state changes to be excluded
 
 twbysleep=struct(); %initalize structure
 
 % Sepearte twitches by each sleep_yn column. Puts all the twitches in one
 % sleep_yn column together
 for s = 1:sleep_num 
     currow=['s' num2str(s)];
     twbysleep.(currow) = struct();
     t_row =[];
     for t = 1:length(twitch.onset)
         tord=twitch.sleepyn_ordinal(t);
         tw=twitch.onset(t);
     if tord == s
         t_row = [t_row tw];
     else
     end
     end 
     twbysleep.(currow) = t_row;
 end
 
 ITIsleep =[]
 ITIsleep_t =[]
 
 %calculate ITIs by each sleep yn column
 for s = 1:sleep_num 
     sleep_yn = sleep_state.yn(s);
     currow=['s' num2str(s)];
     twmatrix = twbysleep.(currow);
     ITI_matrix=[];
     time_matrix=[];
     sleepord.(currow).ITI = struct();
     sleepord.(currow).time = struct();
     if (length(twmatrix))> 0
        for t = 1:(length(twmatrix)-1)
            t1 = twmatrix (t);
            t2 = twmatrix (t+1);
            ITI_matrix(t) = t2-t1;
            time_matrix(t) = (t1+t2)/2;
        end
     else
     end
     sleepord.(currow).ITI = ITI_matrix;
     sleepord.(currow).time = time_matrix;
     
     if sleep_yn == 1
         ITIsleep = [ITIsleep ITI_matrix];
         ITIsleep_t = [ITIsleep_t time_matrix];
     else
     end
 end
 
 ITI= struct()
 ITI.Sleep=struct()
 ITI.Sleep.ITI = ITIsleep';
 ITI.Sleep.Time = ITIsleep_t';
 
 %% Seperate twitch by state ordinal and then calculate ITIs
 %This allows ITIs to be calculated for each state seperately
 twbyordinal=struct();  
 for o = 1:num 
     currow=['O' num2str(o)];
     twbyordinal.(currow) = struct();
     t_row =[];
     for t = 1:length(twitch.onset)
         tord=twitch.ordinal(t);
         tw=twitch.onset(t);
     if tord == o
         t_row = [t_row tw];
     else
     end
     end 
     twbyordinal.(currow) = t_row;
 end

%Initialize
ITIREM = []; 
ITINREM = []; 
ITIN1 = []; 
ITIN2 = []; 
ITIN3 = []; 
ITIREM_t =[]; 
ITINREM_t = []; 
ITIN1_t = []; 
ITIN2_t = []; 
ITIN3_t = []; 

%Calculate ITI in each ordinal sleep subset
for o = 1:num 
     sleep = discrete_state.state(o);
     currow=['O' num2str(o)];
     ITI.ordinal.(currow).ITI = struct();
     ITI.ordinal.(currow).time = struct();
     twmatrix = twbyordinal.(currow);
     ITI_matrix=[];
     time_matrix=[];
     if (length(twmatrix))> 0
        for t = 1:(length(twmatrix)-1)
            t1 = twmatrix (t);
            t2 = twmatrix (t+1);
            ITI_matrix(t) = t2-t1;
            time_matrix(t) = (t1+t2)/2;
        end
     else
     end
     ITI.ordinal.(currow).ITI = ITI_matrix;
     ITI.ordinal.(currow).time = time_matrix;
     
     if sleep == 0
         ITIREM = [ITIREM ITI_matrix];
         ITIREM_t = [ITIREM_t time_matrix];
     else
     end
     
     if sleep == -1
         ITIN1 = [ITIN1 ITI_matrix];
         ITIN1_t = [ITIN1_t time_matrix];
     else
     end
     
     if sleep == -2
         ITIN2 = [ITIN2 ITI_matrix];
         ITIN2_t = [ITIN2_t time_matrix];
     else
     end
     
     if sleep == -3
         ITIN3 = [ITIN3 ITI_matrix];
         ITIN3_t = [ITIN3_t time_matrix];
     else
     end
     
     if sleep == -1 | sleep == -2 | sleep == -3
        ITINREM = [ITINREM ITI_matrix];
        ITINREM_t = [ITINREM_t time_matrix];
     else
     end 
end
 
%Save data in structures
 ITI.bystate.REM.ITI = ITIREM'; 
 ITI.bystate.NREM.ITI = ITINREM'; 
 ITI.bystate.N1.ITI = ITIN1'; 
 ITI.bystate.N2.ITI = ITIN2'; 
 ITI.bystate.N3.ITI = ITIN3'; 
 
 ITI.bystate.REM.Time = ITIREM_t'; 
 ITI.bystate.NREM.Time = ITINREM_t'; 
 ITI.bystate.N1.Time = ITIN1_t'; 
 ITI.bystate.N2.Time = ITIN2_t'; 
 ITI.bystate.N3.Time = ITIN3_t'; 
 
 %% Calculate the Total Time in each state

%Initalize Variables
discrete_state.bystate_dur.Sleep = [];
discrete_state.bystate_dur.REM=[];
discrete_state.bystate_dur.N1=[];
discrete_state.bystate_dur.N2=[];
discrete_state.bystate_dur.N3=[];
 
%Loops through each state row 
for b = 1:num
     sleep = discrete_state.state(b, :);
     dur = discrete_state.dur(b);
     
     
     if sleep == 0 | sleep == -1| sleep == -2 | sleep == -3 %if twitch is during any sleep
         discrete_state.bystate_dur.Sleep = [discrete_state.bystate_dur.Sleep dur];
     else 
     end
     
     if sleep == 0
         discrete_state.bystate_dur.REM = [discrete_state.bystate_dur.REM dur];
     else
     end
     
     if sleep == -1
         discrete_state.bystate_dur.N1 = [discrete_state.bystate_dur.N1 dur];
     else
     end
     
     if sleep == -2
         discrete_state.bystate_dur.N2 = [discrete_state.bystate_dur.N2 dur];
     else
     end
     
     if sleep == -3
         discrete_state.bystate_dur.N3 = [discrete_state.bystate_dur.N3 dur];
     else
     end
end

%Sum of all of the bouts of each state to find total duration in each state
Summary.StateDurations.REM = sum(discrete_state.bystate_dur.REM);
Summary.StateDurations.N1 = sum(discrete_state.bystate_dur.N1);
Summary.StateDurations.N2  = sum(discrete_state.bystate_dur.N2);
Summary.StateDurations.N3  = sum(discrete_state.bystate_dur.N3);
Summary.StateDurations.NREM = Summary.StateDurations.N1 + Summary.StateDurations.N2 + Summary.StateDurations.N3;
Summary.StateDurations.Sleep = Summary.StateDurations.REM + Summary.StateDurations.NREM;

%% Calculates Summary Twitch Variables in each state

%Calc Average Twitch Rate #of twitches/min
Summary.TwitchRate.Sleep = (Summary.TwitchTotal.Sleep/Summary.StateDurations.Sleep)*60;
%Calc Average Twitch Duration
Summary.AverageTwitchDur.Sleep = mean(twitch.dur);

%REM
if Summary.StateDurations.REM >0
Summary.TwitchRate.REM = (Summary.TwitchTotal.REM/Summary.StateDurations.REM)*60;
Summary.AverageTwitchDur.REM = mean(twitch.bystate.REM.dur);
else
    Summary.TwitchRate.REM = NaN;
    Summary.AverageTwitchDur.REM = NaN;
end

%N1
if Summary.StateDurations.N1 >0
Summary.TwitchRate.N1 = (Summary.TwitchTotal.N1/Summary.StateDurations.N1)*60;
Summary.AverageTwitchDur.N1 = mean(twitch.bystate.N1.dur);
else
    Summary.TwitchRate.N1 = NaN;
    Summary.AverageTwitchDur.N1 = NaN;
end

%N2
if Summary.StateDurations.N2 >0
Summary.TwitchRate.N2 = (Summary.TwitchTotal.N2/Summary.StateDurations.N2)*60;
Summary.AverageTwitchDur.N2 = mean(twitch.bystate.N2.dur);
else
    Summary.TwitchRate.N2 = NaN;
    Summary.AverageTwitchDur.N2 = NaN;
end

%N3
if Summary.StateDurations.N3 >0
Summary.TwitchRate.N3 = (Summary.TwitchTotal.N3/Summary.StateDurations.N3)*60;
Summary.AverageTwitchDur.N3 = mean(twitch.bystate.N3.dur);
else
    Summary.TwitchRate.N3 = NaN;
    Summary.AverageTwitchDur.N3 = NaN;
end

%NREM
if Summary.StateDurations.NREM >0
Summary.TwitchRate.NREM = (Summary.TwitchTotal.NREM/Summary.StateDurations.NREM)*60;
Summary.AverageTwitchDur.NREM = mean(twitch.bystate.NREM.dur);
else
    Summary.TwitchRate.NREM = NaN;
    Summary.AverageTwitchDur.NREM = NaN;
end

%% Calculates the Total number of REMs
Summary.TotalREMs = length(rem.onset);
end

