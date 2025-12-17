function Rand = expected_Shuffle(startTimes, endTimes, nShuffles, jitterValue)

% Description:

% Creates a N x M matrix (Rand) of start times, where:
% N = nShuffles 
% M = length(startTimes)
%
% Inputs:
% startTimes = vector of onset times (s) of the event you want to shuffle
% endTimes = vector of offset times (s) of the event you want to shuffle
% nShuffles = How many shuffles of randomized data you want (usually 1000)
% jitterValue = How much do you want to jitter the randomized data by
%               (usually set this to median ITI or ISI)
%
% Used in Christiansen et al., 2026

if size(startTimes,1) == 1
    startTimes = startTimes';
end

if size(endTimes,1) == 1
    endTimes = endTimes';
end

Dur=endTimes-startTimes;
durm1=Dur(2:(length(Dur)));

FirstEvent = startTimes(1);
FirstDur= Dur(1);
AllIntervals = diff(startTimes);

nAllIntervals = length(AllIntervals)+1;

RandStart = zeros(nShuffles, nAllIntervals);
RandDur = zeros(nShuffles, nAllIntervals);
RandMid= zeros(nShuffles, nAllIntervals);

for iRandPerm = 1:nShuffles
    randper=randperm(nAllIntervals-1);
    RandITIs = AllIntervals(randper); %randomize the ITIs
    randdurm1=durm1(randper);
    Rand = [FirstEvent; RandITIs]; %first event plus randomized ITIs
    Rand_dur = [FirstDur; randdurm1];
    cRand = cumsum(Rand); %cumulative sum starts at first event and adds the randomized ITIs iteratively
    adjit = rand*(jitterValue*2); 
    subjit = rand*(jitterValue*2); 
    cRand_jit = cRand + adjit - subjit; %adds or subtracts a randomized jitter
    RandStart(iRandPerm,:) = cRand_jit;
    RandDur(iRandPerm, :) = Rand_dur; 
    midpoint = cRand_jit + (Rand_dur/2);
    RandMid(iRandPerm, :)=midpoint;
end

Rand=struct();
Rand.Start=RandStart;
Rand.Dur=RandDur;
Rand.Mid=RandMid;

end %end function

