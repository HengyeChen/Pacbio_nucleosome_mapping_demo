function NDR_num_len_v2(name,r_index)
%% initiate parameters
%fname = '/storage/home/h/hxc62/work/';%define the home folder
%addpath(fname);
%fname = 'E:/d/Dropbox/';%set the home folder
%list of different pacbio seq samples
%list = [{'WTall'};{'IR2C'};{'NR3Cn'};{'OR3Cn'};{'TRall'};{'NR3C'};{'OR3C'};{'TR'};{'17R'};{'WR2C'};{'N1R3'};{'N3R3'};{'O2R3'};{'O3R3'}];

% if nargin == 1
    border = 73;% set the distance between the boundary and the dyad
    linker = 20;% set the length of standard linker. When predicted linker is short than this value, set a penalty for it. The minmal length of linker is 5bp shorter than this value.
% elseif nargin == 2
%     linker = 20;
% end

motifs = load('ref\motif_pos_v2.mat');%read motif positions
roi = 250:803;%region of interest on each sequence
%% calcualte and save ndr length information
    %name = list{k,1};% the name of sample
    k=1;
    data = load(['nuc_prediction\tt\nuc_pred_TTxx_',name,'.mat']);%load referred nucleosome positioning info
    tic;
    rate = zeros(1,169);%
    thub= [0,0];% initiate threshold value
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gcn = [19,17];% gc number in the two background sequences
    %% calculate methylation threshold based on the methylation level in background sequences
    for i = 144:145
        met_all = data.data(i).met;% read methylation matrix
        mm = met_all(:,roi);% methylation in the region of interest
        mm(mm==-1)=0;% convert GT to 0, GC is 1 
        xxx = sum(mm,2)/gcn(1,i-143);% calculate the average methylation level in each read
        ds = fitdist(xxx,'Normal');% fit the methylation level in each read to normal distribution
        %thub(1,i-143) = ds.mu + 1.645*ds.sigma; % set the threshold at 95%
        thub(1,i-143) = ds.mu + 1.282*ds.sigma;% set the threshold at 90%
    end
    thu = mean(thub);% the final threshold is the average of the two background sequence
    %% calculate NDR length of each sequence 
    for ri = 1:length(r_index)
        i = r_index(ri);
        %set the position of motifs
        if i ~= 144 && i~= 145
            mot_pos = motifs.data(i).pos;
        else
            mot_pos = [488;508];
        end
        %read inferred nucleosome positioning data
        peak_all = data.data(i).pred;
            met_all = data.data(i).met;
        [len,~] = size(peak_all);
        %%%%%%%%%%%%%%%%%%%%%%%
        motu = mot_pos(1,1);%the upstream boundary of motifs
        motd = mot_pos(end,end);%the downstream boundary of motifs
        %% pre-define the size of matrix for NDR info storage
        ndr = zeros(len,2);
        ndr_more = zeros(len,1);
        symm = zeros(len,2);
        ndr_mm = zeros(len,1);
        
        %% calculate NDR length in each read
        for j = 1:len
            metc = met_all(j,:);
            peak = peak_all(j,:);
            n=0;
            m=0;
            if sum(peak(1,motu:motd)<=0)>0%when motif region is not occupied by nucleosome
                ndr(j,1) = 1;
                ndr_more(j,1) = 1;
                ndr_mm(j,1) = 1;
                if peak(1,motu)>0 && peak(1,motd)>0% when upstream and downstream of motif region are occupied by nucleosome
                    ndr(j,2) = sum(peak(1,motu:motd)<=0);
                elseif peak(1,motu)>0 && peak(1,motd)<=0% when upstream of motif region is occupied, but downstream is open
                    for n = 1:1200
                        down = motd+n;
                        if peak(1,down) > 0 || down == 1278
                            break
                        else
                        end
                    end
                    ndr(j,2) = sum(peak(1,motu:motd)<=0)+n-1;
                elseif peak(1,motu)<=0 && peak(1,motd)>0% when downstream of motif region is occupied, but upstream is open
                    for m = 1:1200
                        up = motu-m;
                        if peak(1,up) > 0 || up == 1
                            break
                        else
                        end
                    end
                    ndr(j,2) = sum(peak(1,motu:motd)<=0)+m-1;
                else% when entire motif region is nucleosome free
                    for n = 1:1200 % find how long is the NDR in the downstream of motif
                        down = motd+n;
                        if peak(1,down) > 0 || down == 1278
                            break
                        else
                        end
                    end
                    ndr(j,2) = sum(peak(1,motu:motd)<=0)+n-1;
                    for m = 1:1200
                        up = motu-m;
                        if peak(1,up) > 0 || up == 1
                            break
                        else
                        end
                    end
                    ndr(j,2) = sum(peak(1,motu:motd)<=0)+n+m-2;         
                end 
                symm(j,:) = [m n];
            elseif sum(met_all(j,541:574)) >1
                ndr_more(j,1) = 1;
            else
                    ndr(j,1) = 0;
                    ndr(j,2) = 0;
                    ndr_more(j,1) = 0;
                    symm(j,:) = [0 0]; 
            end
            if ndr(j,1) == 0 && (sum(metc(1,roi)==1)/sum(metc(1,roi)~=0))>thu
                
                ndr_mm(j,1) = 1;
            end
        end
        %% save the data in the structure 'all'
        all(k).info(i).num = len;
        all(k).info(i).ndr_rate = sum(ndr(:,1))/len;
        rate(1,i) = sum(ndr(:,1))/len;
        ndrs = ndr(:,2);
        ndrs_only = sort(ndrs(ndrs~=0),'descend');% sort the NDR length
        mndrso = mean(ndrs_only((1+round(len/20)):end,1));%calculate average NDR length of reads containing NDRs(top %5 longest reads are removed because when NDRs are short, this small population of overmethylated reads will significantly affect the average length)
        all(k).info(i).ndr_mean = sum(ndrs)/len;% average NDR length of all reads
        all(k).info(i).ndr_median = median(ndrs);% median NDR length of all reads
        all(k).info(i).ndrs = ndrs;% NDR length of each read
        all(k).info(i).ndr_only_mean = mndrso;% average NDR length of reads containing NDRs(top %5 reads are removed), this is the average NDR length used in paper
        %ndr_info(i,:) = [i len mndrso sum(ndr(:,1))/len];
    end   
    %th = 140;
    %tmot(k).lp = NDR_mean_std_calculation_3motif(ndr_info,th,fname);
    %tmot(k).lp = NDR_mean_std_calculation_comb(ndr_info,th);
save_data(['nuc_prediction\tt\pred_NDR_lp_',name,'.mat'],all);
end
