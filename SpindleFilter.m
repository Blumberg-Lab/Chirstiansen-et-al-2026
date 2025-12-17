function [fx,fTimes] = SpindleFilter(x,Times,Fs)
% [fx,fTimes] = SpindleFilter(x,Times,Fs)
%
% Filters LFP data x in spindle frequencies and removes phase shift
% introduced by the filter.
%
% Input: x      Waveform Values
%        Times  Waveform Times
%        Fs     Waveform sampling rate (Hz)
%
% Output:fx     Waveform values after filtering
%        fTimes Adjusted Waveform times (Eliminating phase shift).
%
% Anotated 6/15/2021. 
%
% Written 9/23/2020 by Jimmy Dooley (james-c-dooley@uiowa.edu)
%
% CC-by-4.0 to Jimmy Dooley
% Used in Sokoloff et al., 2021.
%
% Cite as:
% Sokoloff G, Dooley JC, Glanz RM, Yen RY, Hickerson MM, Evans L, Laughlin
%    HM, Apfelbaum KS, and Blumberg MS (2021). Twitches emerge postnatally 
%    during quiet sleep in human infants and are synchronized with sleep 
%    spindles. Current Biology. https://doi.org/10.1016/j.cub.2021.05.038

% TC updated 5/21/24: stopband changed to 10-16 Hz and passband 11-15 Hz 
% in accordance with McClain et al., 2016 & D'Atri et al., 2018 

persistent Hd;

if isempty(Hd)
    
    Fstop1 = 10;        % First Stopband Frequency
    Fpass1 = 11;        % First Passband Frequency
    Fpass2 = 15;        % Second Passband Frequency
    Fstop2 = 16;        % Second Stopband Frequency
    Astop1 = 60;        % First Stopband Attenuation (dB)
    Apass  = 1;         % Passband Ripple (dB)
    Astop2 = 60;        % Second Stopband Attenuation (dB)
    
    h = fdesign.bandpass('fst1,fp1,fp2,fst2,ast1,ap,ast2', Fstop1, Fpass1, ...
        Fpass2, Fstop2, Astop1, Apass, Astop2, Fs);
    
    Hd = design(h, 'equiripple', ...
        'MinOrder', 'any');
    
    
    set(Hd,'PersistentMemory',true);
    
end

delay = round(length(Hd.States)/2); % Phase shift is equal to half the order of the filter

y = filter(Hd,x);
fx = y(delay+1:end);
fTimes = Times(1:end-delay);





