function all = NDR_num_len_v2(predicted,motif,is_bkg,bkg_pseudo_pos,border,length_proportion)
%% this function calculate the length and proportion of nucleosome-depleted regions.
%predicted is the output of nuc_predict_v2.m. It contains the nucleosome prediction matrix.
%motif is the motif information file which is a mat file containing the position, sequence, and name of each TF motif
%is_bkg determines whether the seqence is a background sequence. 1 is yes. 0 is no.
%bkg_pseudo_pos is the pseudo position of the background sequence
%border is the nucleosome protection length on one side of the nucleosome dyad
%length_proportion is the output file containing NDR length and proportion

tic
%% calculate NDR length of each sequence
%set the position of motifs
if is_bkg == 1
    mot_pos = bkg_pseudo_pos;
else
    mot_pos = load(motif).mot.pos;
end
%read inferred nucleosome positioning data
peak_ndf = predicted.pred;
met_ndf = predicted.met;
[len,~] = size(peak_ndf);
%%%%%%%%%%%%%%%%%%%%%%%
motu = mot_pos(1,1);%the upstream boundary of motifs
motd = mot_pos(end,end);%the downstream boundary of motifs
%% pre-define the size of matrix for NDR info storage
ndr = zeros(len,2);
%% calculate NDR length in each read
for j = 1:len
    metc = met_ndf(j,:);
    peak = peak_ndf(j,:);
    if sum(peak(1,motu:motd)<=0)>0%when motif region is not occupied by nucleosome
        ndr(j,1) = 1;
        if peak(1,motu)>0 && peak(1,motd)>0% when upstream and downstream of motif region are occupied by nucleosome
            ndr(j,2) = sum(peak(1,motu:motd)<=0);
        elseif peak(1,motu)>0 && peak(1,motd)<=0% when upstream of motif region is occupied, but downstream is open
            for n = 1:(length(metc)-border)
                down = motd+n;
                if peak(1,down) > 0 || down == length(metc)
                    break
                else
                end
            end
            ndr(j,2) = sum(peak(1,motu:motd)<=0)+n-1;
        elseif peak(1,motu)<=0 && peak(1,motd)>0% when downstream of motif region is occupied, but upstream is open
            for m = 1:(length(metc)-border)
                up = motu-m;
                if peak(1,up) > 0 || up == 1
                    break
                else
                end
            end
            ndr(j,2) = sum(peak(1,motu:motd)<=0)+m-1;
        else% when entire motif region is nucleosome free
            for n = 1:(length(metc)-border) % find how long is the NDR in the downstream of motif
                down = motd+n;
                if peak(1,down) > 0 || down == length(metc)
                    break
                else
                end
            end
            ndr(j,2) = sum(peak(1,motu:motd)<=0)+n-1;
            for m = 1:(length(metc)-border)
                up = motu-m;
                if peak(1,up) > 0 || up == 1
                    break
                else
                end
            end
            ndr(j,2) = sum(peak(1,motu:motd)<=0)+n+m-2;
        end
    else
        ndr(j,1) = 0;
        ndr(j,2) = 0;
    end
end
%% save the data in the structure 'all'
all.info(1).num = len;
all.info(1).ndr_rate = sum(ndr(:,1))/len;
ndrs = ndr(:,2);
ndrs_only = sort(ndrs(ndrs~=0),'descend');% sort the NDR length
mndrso = mean(ndrs_only((1+round(len/20)):end,1));%calculate average NDR length of reads containing NDRs(top %5 longest reads are removed because when NDRs are short, this small population of overmethylated reads will significantly affect the average length)
all.info(1).ndr_mean = sum(ndrs)/len;% average NDR length of all reads
all.info(1).ndr_median = median(ndrs);% median NDR length of all reads
all.info(1).ndrs = ndrs;% NDR length of each read
all.info(1).ndr_only_mean = mndrso;% average NDR length of reads containing NDRs(top %5 reads are removed), this is the average NDR length used in paper
%ndr_info(i,:) = [i len mndrso sum(ndr(:,1))/len];
save_data(length_proportion,all);
disp('NDR length and proportion are calculated');
toc
end

