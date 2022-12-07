function alignment = align_met_all_v2(seqs_c,refs,motifs,variable_region,align_accuracy,output)
%%This function align reads to the reference sequence
%seqs_c is the output of rev2fw_header_v2.m. It should contains raw reads and reverse complimentary reads
%refs is the reference file which is a txt file.
%motifs is the motif information file which is a mat file containing the position, sequence, and name of each TF motif
%variable_region is the position of the variable region in the reference sequence
%output is the output file
tic;
seqs = seqs_c(:,1);
headers = seqs_c(:,2);
mot_pos = load(motifs);% read positions and sequences of motifs
references = readcell(refs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for r = 1
    ref = references{r,1};%read a reference sequence
    ref_roi = ref(variable_region);% the variable region (the region that contains motifs) in the reference sequence
    count = length(seqs);
    seqs = seqs(~cellfun('isempty',seqs));
    mheader = cell(count,1);
    aligned = char(zeros(count,length(ref)));
    ac = 0;
    for i = 1:count% align each seq to the selected ref
        header = headers{i,1};
        try
            %% align each read to the variable region in the reference
            [~, roi_local] = swalign(ref_roi,seqs{i,1},...%locally align unmatched seqs to the ref seq using Smith-Waterman algorithm
                'Alphabet', 'NT',...%set seq type as nucleotide('NT')
                'ScoringMatrix','NUC44');%use matrix NUC44 to align sequence, only the well-matched part will be aligned
            % If the read can be aligned to the variable region, then align
            % it to the full-length reference sequence. This step can
            % accelerate the alignment
            if sum(roi_local(2,:)== '|' | roi_local(2,:)==':')>80
                %% align each read to the reference sequence
                [~, local_align] = swalign(ref,seqs{i,1},...%locally align unmatched seqs to the ref seq using Smith-Waterman algorithm
                    'Alphabet', 'NT',...%set seq type as nucleotide('NT')
                    'ScoringMatrix','NUC44');%use matrix NUC44 to align sequence, only the well-matched part will be aligned
                %% In this module the indel mutations in the sequencing reads are removed, and the output data aseq has the same length as the ref seq
                local_match = strrep(local_align(3,:),'-','N');%replace all the gaps in the alignment to N
                [~,align] = nwalign(ref,local_match,...%globally align realigned seqs to the ref using Needleman-Wunsch algorithm
                    'Alphabet', 'NT',...%set seq type as nucleotide('NT')
                    'ScoringMatrix','NUC44');%use matrix NUC44 to align sequence, the whole sequence will be aligned
                
                alignr = [align(1,:);align(3,:)];%get the aligned sequences from the alignment output data
                [~,w] = size(alignr);% get the length of the aligned sequence
                for width = 1:w
                    if alignr(2,width)=='-'%replace the gaps in the sequence by N
                        alignr(2,width)=strrep(alignr(2,width),'-','N');
                    end
                    if alignr(1,width)=='-'%replace the insert mutations in the sequence by x
                        alignr(2,width)='x';
                    end
                end
                
                aseq=regexprep(alignr(2,:),'x','');% delete all x in the sequence, so that the length of aseq is same as the reference seq
                %% determine whether the sequence is well aligned or not.
                [~, align2] = nwalign(ref,aseq,...% globally align the trimmed aseq to ref to get the similarity
                    'Alphabet', 'NT',...
                    'ScoringMatrix','NUC44');
                
                sim_all = sum(align2(2,:)== '|' | align2(2,:)==':');%count how many nt are matched in the whole seq
                sim_v = sum(align2(2,variable_region)== '|' | align2(2,variable_region)==':');%count how many nt are matched in the variable region
                %if the overall similarity is above 90% and similarity at variable region is higher than 99%, add this sequence to 'aligned'
                if sim_all > align_accuracy*length(ref) && sim_v > 0.99*length(variable_region)
                    real = motif_align(aseq,mot_pos,r);
                    if real == 1
                        ac = ac+1;
                        aligned(ac,:) = aseq;
                        mheader{ac,1} = header;
                        seqs{i,1}=[];% Aligned reads are removed from the list, so it will not be aligned to next reference sequence.
                    end
                end
            end
        catch
        end
    end
    if ac~=0
        alignment.reads = aligned(1:ac,:); %input the aligned seq to the struct, r is the index of the ref
        alignment.header = mheader(~cellfun('isempty',mheader));%remove empty variables in the cell array mheader
    end
    disp(r);
end
save(output, 'alignment');%save all rematch seqs
disp('alignment completed');
toc;
end