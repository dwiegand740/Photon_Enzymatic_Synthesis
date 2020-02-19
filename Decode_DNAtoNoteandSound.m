%% Decoding (From DNA to Ternary,MIDI note and play)
%Howon Lee - Harvard Medical School, Church lab, howon251@gmail.com,
clear all
clc

%% Sequencing result
Seqresult = [
    'GATATGCT'
    'GATGCTAG'
    'GTCTGCAT'
    'GAGACACT'
    'GAGCATGA'
    'GACTATAG'
    'GACGCATC'
    'GACACGAG'
    'GTCGATAG'
    'GTCATGCA'
    'GAGTATAG'
    'GATCGCGA'];
%% DNA to Ternary conversion
for se = 1:length(Seqresult)
    tpseq = Seqresult(se,:);
    m = 1;
    for tp = 1:length(tpseq)-1
        %from G
        
        if tpseq(tp) == 'G' & tpseq(tp+1) == 'A'
            dt(m) = '0';
            nuclnote(m,:) = ['G','A'];
            m = m+1;
            
        elseif tpseq(tp) == 'G' & tpseq(tp+1) == 'C'
            dt(m) = '2';
            nuclnote(m,:) = ['G','C'];
            m = m+1;
            
        elseif tpseq(tp) == 'G' & tpseq(tp+1) == 'T'
            dt(m) = '1';
            nuclnote(m,:) = ['G','T'];
            m = m+1;
            
            %from T
        elseif tpseq(tp) == 'T' & tpseq(tp+1) == 'A'
            dt(m) = '1';
            nuclnote(m,:) = ['T','A'];
            m = m+1;
            
        elseif tpseq(tp) == 'T' & tpseq(tp+1) == 'C'
            dt(m) = '0';
            nuclnote(m,:) = ['T','C'];
            m = m+1;
            
        elseif tpseq(tp) == 'T' & tpseq(tp+1) == 'G'
            dt(m) = '2';
            nuclnote(m,:) = ['T','G'];
            m = m+1;
            
            %from A
        elseif tpseq(tp) == 'A' & tpseq(tp+1) == 'C'
            dt(m) = '1';
            nuclnote(m,:) = ['A','C'];
            m = m+1;
            
        elseif tpseq(tp) == 'A' & tpseq(tp+1) == 'G'
            dt(m) = '0';
            nuclnote(m,:) = ['A','G'];
            m = m+1;
            
        elseif tpseq(tp) == 'A' & tpseq(tp+1) == 'T'
            dt(m) = '2';
            nuclnote(m,:) = ['A','T'];
            m = m+1;
            
            %from C
        elseif tpseq(tp) == 'C' & tpseq(tp+1) == 'A'
            dt(m) = '2';
            nuclnote(m,:) = ['C','A'];
            m = m+1;
            
        elseif tpseq(tp) == 'C' & tpseq(tp+1) == 'G'
            dt(m) = '1';
            nuclnote(m,:) = ['C','G'];
            m = m+1;
            
        elseif tpseq(tp) == 'C' & tpseq(tp+1) == 'T'
            dt(m) = '0';
            nuclnote(m,:) = ['C','T'];
            m = m+1;
            
        end
    end
    
    infoTer(se,:) = dt;% ternary
    
    
end % se

%% Ternary to Decimal conversion
infoTerInd = infoTer(:,1:3);
infoTerNote = infoTer(:,4:6);
infoTerTem = infoTer(:,end);

infoDecInd = base2dec(infoTerInd,3);
infoDecNote = base2dec(infoTerNote,3);
infoDecTem = base2dec(infoTerTem,3);

% sort by index
infoDec = [infoDecInd infoDecNote infoDecTem];

for indso = 1:length(infoDec)
    infoDecInd(indso)
    Note(infoDecInd(indso)+1) = infoDecNote(indso);
    Tempo(infoDecInd(indso)+1) = infoDecTem(indso);
        
    
end% for indso


%% Midi play and wav file generation

%constant
dur = 0.65; %second
oct = 12*5; %default octave

note = Note;
durationT = Tempo;
wx = [];

for n = 1:length(note)
    
    if note(n)>25 % Mute
        note(n) = -oct+20;
    else
        note(n) = note(n);
    end
    
    midinum = oct+note(n); % Recover original Midi number
    dura = dur/4*(durationT(n)+1); % Play duration
    freq = 8.1757989156*2^(midinum/12);% Frequency of sound
    
    fs = 44100;% Sampling frequency
    t = 0:1/fs:dura;
    x = sin(2*pi*freq*t);
    soundsc(x, fs);% Sound generation
    wx = [wx x];
    pause(dur/4*(durationT(n)+1)) 
    
    
end


filename = 'SuperMario_Howon.wav';
audiowrite(filename,wx,fs);


