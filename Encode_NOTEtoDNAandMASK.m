%% Super Mario from note number to DNA conversion
%Howon Lee - Harvard Medical School, Church lab, howon251@gmail.com,
%617-817-5122, 12/22/2019
clear all
clc
%constant
dur = 0.65;%second
oct = 12*5;
% use only octave 5,6,7 (note number 60~85, C5~C#7, note number 86 means MUTE)
% use only 1/4, 2/4/, 3/4 note duration indicated by 0 (quarter note), 1 (half note), 2 each

%% Input note information
noteO = [76,76,86,76,86,72,76,86,79,86,67,86];% input note number
note = noteO-oct; % substract default octave
note_input = note;
duration = [0,0,0,0,0,0,0,0,0,2,0,2];% input note duration

%% Conversion

%Original data
ind = [0:length(noteO)-1];
mes = note_input;
duration = duration;
%Decimal to Ternary conversion
tr_ind = dec2base(ind,3);
tr_mes = dec2base(mes,3);
tr_duration = dec2base(duration,3);

data = [tr_ind tr_mes tr_duration];

%% Ternary to DNA conversion 
for n=1:size(data,1)
    
    indata = data(n,:);
    %initial sequence "G"
    seq = 'G';
    for m=1:length(indata)
        %from G
        if seq(end)=='G' & indata(m)=='0'
            seq = [seq,'A'];
        elseif seq(end)=='G' & indata(m)=='1'
            seq = [seq,'T'];
        elseif seq(end)=='G' & indata(m)=='2'
            seq = [seq,'C'];
            %from T
        elseif seq(end)=='T' & indata(m)=='0'
            seq = [seq,'C'];
        elseif seq(end)=='T' & indata(m)=='1'
            seq = [seq,'A'];
        elseif seq(end)=='T' & indata(m)=='2'
            seq = [seq,'G'];
            %from A
        elseif seq(end)=='A' & indata(m)=='0'
            seq = [seq,'G'];
        elseif seq(end)=='A' & indata(m)=='1'
            seq = [seq,'C'];
        elseif seq(end)=='A' & indata(m)=='2'
            seq = [seq,'T'];
            %from C
        elseif seq(end)=='C' & indata(m)=='0'
            seq = [seq,'T'];
        elseif seq(end)=='C' & indata(m)=='1'
            seq = [seq,'G'];
        elseif seq(end)=='C' & indata(m)=='2'
            seq = [seq,'A'];
        else
            seq = seq;
        end
        
    end
    
    template(n,:) = seq;
    
end


%% Mask generation 

%base mask (single spot) generation

x_pix = 1080;
y_pix = 1920;

mask = zeros(x_pix,y_pix,12);

%position, radius
posX = x_pix/2;
posY = y_pix/2;
radi = 100/2;
fill = radi*8;

x = [posX-fill, posX, posX+fill];
y = [posY-fill-fill/2, posY-fill/2,posY+fill/2, posY+fill+fill/2];

for xn=1:length(x)
    for yn=1:length(y)
        for n=1:size(mask,1)
            for m =1:size(mask,2)
                distance = sqrt((n-x(xn))^2+(m-y(yn))^2);
                if distance<radi
                    mask(n,m,(xn-1)*4+yn)=1;
                else
                    mask(n,m,(xn-1)*4+yn)=mask(n,m,(xn-1)*4+yn);
                end
            end
        end
    end
end

% Synthesis mask generation using base mask
for tn = 2:size(template,2)
    
    % A, T, G, C synthesis of tn cycle
    
    % A mask
    tmseq = template(:,tn);
    mask_A = zeros(1080,1920); %DMD resolution
    for amask =1:length(tmseq)
        if tmseq(amask)=='A'
            mask_A = mask_A + mask(:,:,amask);
        else
            mask_A = mask_A;
        end
    end
    filename = strcat('0',num2str(tn-1),'_','0A','_','mask.bmp');
    imwrite(mask_A,filename)
    
    % T mask
    mask_T = zeros(1080,1920);
    for tmask =1:length(tmseq)
        if tmseq(tmask)=='T'
            mask_T = mask_T + mask(:,:,tmask);
        else
            mask_T = mask_T;
        end
    end
    filename = strcat('0',num2str(tn-1),'_','1T','_','mask.bmp');
    imwrite(mask_T,filename)
    
    % G mask
    mask_G = zeros(1080,1920);
    for gmask =1:length(tmseq)
        if tmseq(gmask)=='G'
            mask_G = mask_G + mask(:,:,gmask);
        else
            mask_G = mask_G;
        end
    end
    filename = strcat('0',num2str(tn-1),'_','2G','_','mask.bmp');
    imwrite(mask_G,filename)
    
    % C mask
    mask_C = zeros(1080,1920);
    for cmask =1:length(tmseq)
        if tmseq(cmask)=='C'
            mask_C = mask_C + mask(:,:,cmask);
        else
            mask_C = mask_C;
        end
    end
    filename = strcat('0',num2str(tn-1),'_','3C','_','mask.bmp');
    imwrite(mask_C,filename)
    
    
end

save Midi_NOTEtoDNAandMASK
