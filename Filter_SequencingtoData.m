%% "Photon-directed Multiplexed Enzymatic DNA Synthesis for Molecular Digital Data Storage"
%% Sequence filtering - Howon Lee (howon251@gmail.com)
% Initial boundary condition
% 1)Sequence of Illumina adaptor
% 2)# of transition
% 3)Sequence of locational barcode
% 4)# of strands synthesized
% 5) Starting nucleotide : G

clear all
clc

%% 1. Loading Fastq
mario31 = fastqread('C:\Users\Wyss User\Desktop\TdT\Sequencing\SuperMario\Howon_Library_4_Read1_10102019.fastq');

%% 2. Prescreening with illumina adaptor sequence (Filter #1)
% st : starting adaptor
% ed : ending adaptor

st = 'TAGTGTGCTTCGGACC';
ed = 'CGACTGAACCCAAGCA';

% Sort out sequences with starting and ending universal adaptor sequences only
% Using Smith-Waterman algorithm

m=1;% # of sorted sequences
trun_data={};% sorted sequences
for n =1:size(mario31,2)
    
    test = mario31(n).Sequence;
    
    % for starting adaptor sequence
    [score,alignment,start] = swalign(test,st);
    
    if score>=30 % sort out sequence with matching score over 30 only
        initial = start(1);
        mi = 1;
    else
        initial = 9999;
        mi = 0;
    end
    
    data_start = initial+length(st);
    
    %for ending adpator sequence
    [scored,alignmentd,startd] = swalign(test,ed);
    if scored>=30
        ending = startd(1);
        me = 1;
    else
        ending = 9999;
        me = 0;
    end
    
    data_end = ending-1;
    
    % Sort out both starting and ending adpator sequences only
    if mi+me==2 & isempty(test(data_start:data_end))==0
        trun_data{m} = test(data_start:data_end);
        m=m+1;
    else
        m=m;
    end
end

% Sequences with both starting and ending adpator sequences
trun_data = trun_data';

%% # of sequencing result(total)    : 216252
%  # of sequences with both adaptor : 100760
%  46.6% of sequences have illumina adaptor at both end

%% 3. Extract seqence transition and its homopolymer length information

transition = {}; % Sequence transition information only
transition_length = {}; % and corresponding homopolymer length information

for k = 1:length(trun_data)
    
    read = trun_data{k};
    
    tran = read(1);
    homo_length = [];
    f2 = 1;
    
    for f=2:length(read)
        
        if read(f-1)==read(f)
            tran = tran;
            f2 = f2+1;
        else
            tran = [tran read(f)]; % only record when sequence transition event occur
            homo_length = [homo_length f2];% record each homopolymer length
            f2 = 1;
        end
    end
    
    homo_length = [homo_length f2];
    
    transition{k} = tran;
    transition_length{k} = homo_length;
    
end

transition = transition';
transition_length = transition_length';

%% Subset generation
%Find sequences with exact number of sequence transition
%Find sequences with exact locational barcode

% Location barcode only (with unknown 4~5nt at the end)
temp_ref(1).seq='GAGAXXXXX';
temp_ref(1).length=length(temp_ref(1).seq);

temp_ref(2).seq='GAGTXXXXX';
temp_ref(2).length=length(temp_ref(2).seq);

temp_ref(3).seq='GAGCXXXXX';
temp_ref(3).length=length(temp_ref(3).seq);

temp_ref(4).seq='GACTXXXXX';
temp_ref(4).length=length(temp_ref(4).seq);

temp_ref(5).seq='GACGXXXX';
temp_ref(5).length=length(temp_ref(5).seq);

temp_ref(6).seq='GACAXXXXX';
temp_ref(6).length=length(temp_ref(6).seq);

temp_ref(7).seq='GATCXXXXX';
temp_ref(7).length=length(temp_ref(7).seq);

temp_ref(8).seq='GATAXXXXX';
temp_ref(8).length=length(temp_ref(8).seq);

temp_ref(9).seq='GATGXXXXX';
temp_ref(9).length=length(temp_ref(9).seq);

temp_ref(10).seq='GTCTXXXXX';
temp_ref(10).length=length(temp_ref(10).seq);

temp_ref(11).seq='GTCGXXXXX';
temp_ref(11).length=length(temp_ref(11).seq);

temp_ref(12).seq='GTCAXXXXX';
temp_ref(12).length=length(temp_ref(12).seq);
%% Boolean selection

for ind = 1:12 % # of synthesized strands
    ind
    trnum = length(temp_ref(ind).seq);
    
    
    
    %% Locational barcode matching (Filter #2)
    filter2=0;
    for ie = 1:length(transition)
        
        if length(transition{ie})<5
            filter2(ie) = 0;
        else
            
            if transition{ie}(1:4)==temp_ref(ind).seq(1:4)
                filter2(ie) = 1; %set value 1 when it's true
            else
                filter2(ie) = 0;
            end
        end
    end
    
    %% # of sequence transition (Filter #3)
    filter3=0;
    for io = 1:length(transition_length)
        
        if length(transition_length{io})==trnum
            filter3(io) = 1; % set value 1 when it's true
        else
            filter3(io) = 0;
        end
        
    end
    
    
    filter23 = floor((filter2+filter3)/2);
    filter23_total(:,ind) = filter23;
    
    % Sequence index finder after filter 2 and 3
    f23ind = [];
    for f23=1:length(filter23)
        
        if filter23(f23) == 1
            f23ind = [f23ind f23];
        else
            f23ind = f23ind;
        end
        
    end
    
    %Sequence logo display
    tmps = transition(f23ind);
    seqlogo(tmps)
    
    f2num(ind) = sum(filter2);
    f23num(ind) = sum(filter23);
    
    
end % for ind

%% Perfect Match (Filter #4)
temp_ref(1).seq='GAGACACTC';
temp_ref(1).length=length(temp_ref(1).seq);

temp_ref(2).seq='GAGTATAGC';
temp_ref(2).length=length(temp_ref(2).seq);

temp_ref(3).seq='GAGCATGAC';
temp_ref(3).length=length(temp_ref(3).seq);

temp_ref(4).seq='GACTATAGC';
temp_ref(4).length=length(temp_ref(4).seq);

temp_ref(5).seq='GACGCATC';
temp_ref(5).length=length(temp_ref(5).seq);

temp_ref(6).seq='GACACGAGC';
temp_ref(6).length=length(temp_ref(6).seq);

temp_ref(7).seq='GATCGCGAC';
temp_ref(7).length=length(temp_ref(7).seq);

temp_ref(8).seq='GATATGCTC';
temp_ref(8).length=length(temp_ref(8).seq);

temp_ref(9).seq='GATGCTAGC';
temp_ref(9).length=length(temp_ref(9).seq);

temp_ref(10).seq='GTCTGCATC';
temp_ref(10).length=length(temp_ref(10).seq);

temp_ref(11).seq='GTCGATAGC';
temp_ref(11).length=length(temp_ref(11).seq);

temp_ref(12).seq='GTCATGCAC';
temp_ref(12).length=length(temp_ref(12).seq);

f4ind_total=[];
for ind = 1:12
    ind
    filter4=0;
    for per = 1:length(transition)
        if filter23_total(per,ind)==1 & transition{per} == temp_ref(ind).seq
            filter4(per) = 1;
        else
            filter4(per)=0;
        end
    end% for f4
    
    f4ind=[];
    for f4=1:length(filter4)
        
        if filter4(f4) == 1
            f4ind = [f4ind f4];
        else
            f4ind = f4ind;
        end
        
    end
    
    f4ind_total = [f4ind_total f4ind]; % Perfect match read position
    
    
end % for ind





