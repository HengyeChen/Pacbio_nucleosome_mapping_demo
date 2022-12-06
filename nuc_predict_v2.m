function infom = nuc_predict_v2(matrix_bkg,matrix_ndf,motif_ndf,bkg_pseudo_pos,border,linker,is_bkg,is_remove_nuc,c_limit,roi,pred_matrix,pred_pileup)%This version predict nuc-6 at the canonical place, when no methylated GC at the upstream of the sequence
%% this function predict nucleosomes from the methylation matrix
%matrix_bkg and matrix_ndf are outputs from pbmatrix_v2.m. They contains methylation matrices for background and ndf sequences
%motif_ndf is the motif information file which is a mat file containing the position, sequence, and name of each TF motif
%bkg_pseudo_pos is the pseudo position of the background sequence
%border is the nucleosome protection length on one side of the nucleosome dyad
%linker is the length of linker
%is_bkg determines whether the seqence is a background sequence. 1 is yes. 0 is no.
%is_remove_nuc determines whether nucleosomes are removed based on methylation levels. 1 is yes. 0 is no.
%c_limit is the maximum optimaziation cycles.
%roi is the position of region of interest
%pred_matrix is the output file containing methylation and prediction matrices
%pred_pileup is the output file containing the bulk methylation level and nucleosome occupancy

tic;
motn = load(motif_ndf);%read motif positions
%%%%%%%%%%%%%%%%%% generate a cumulative distribution curve, and use it as the scoring matrix for a mononucleosome
a = -20; b = border;
x = a:1:b;%set the range of the cumulative distribution curve
m = (a+b)/2;%calculate the mean of a and b, which also is the middle point of x
s = 10; % set the standard deviation of the curve
pd = makedist('Normal','mu',m,'sigma',s);%creates a probability Normal distribution object in which 'mu'(the mean) is m, and 'sigma'(standard deviation) is s
f = cdf(pd,x);%Create a standard normal distribution object with the mean, ? is m and the standard deviation, ? is s. Each value in f correspond to a value in x.
score_matrix = [f(1,22:end), 1, flip(f(1,22:end))];%trim the distribution array, flip the trimmed array, and combine these two arrays. The final array is the score matrix for a nulceosome
lenmat = size(matrix_bkg.C_T_sum_trim,2);%length of a read
bkg_nuc_number = floor(lenmat/(2*b+1+linker));%The maximum number of nucleosomes.
%% nucleosome prediction
dbstop if error
%% set the methylation threshold with the methylation level in background sequences
if is_remove_nuc == 1 
    met_all = matrix_bkg.C_T_sum_trim;
    mm = met_all(:,roi);
    mm(mm==-1)=0;
    xxx = sum(mm,2)/length(find(matrix_bkg.C_T_pos>=roi(1,1) & matrix_bkg.C_T_pos<=roi(end,end)));
    ds = fitdist(xxx,'Normal');
    thu = ds.mu + 1.282*ds.sigma;
end

%% predict nucleosome in each sequence
for i = 1
    mpos = motn.mot.pos;%read the position file
if is_bkg == 1
    mpos = bkg_pseudo_pos;
end
    matrix = matrix_ndf.C_T_sum_trim;%read the trimmed sequence matrix
    [len,~] = size(matrix);%count the number of reads in the matrix
    if len == 0
        continue
    end
    dyad_all = zeros(len,bkg_nuc_number);%create an zero matrix for the positions of dyads of all sequences
    %% predict nucleosome positions
    map = zeros(len,lenmat);%
    nuc_shift_all = zeros(len,lenmat);% create an zero matrix for the shifted nucleosome positions
    score_all_pos = zeros(len,lenmat);
    shift_dyad_all = zeros(len,bkg_nuc_number);% create an zero matrix for the shifted dyad positions
    mot_down = max(max(mpos))+10;%set the boundary of the motifs, downstream + 10bp
    mot_up = min(min(mpos))-10;%set the boundary of the motifs, upstream - 10bp
    pileup_m1 = zeros(len,lenmat);
    gc_pos = matrix_ndf.C_T_pos;
    pileup_met = zeros(3,length(gc_pos));
    seqt = matrix;% create a new array which is same as seq
    seqt(seqt<0) = 0;
    for m = 1:len
        seq = matrix(m,:);%read a single sequence from the matrix
        dyad_f = zeros(1,bkg_nuc_number);%reset dyad positions of all nucleosomes
        score_sum = zeros(1,bkg_nuc_number)-2;% set the default score of each nucleosome as -2
        
        for a = (1+border):(lenmat-border)
            seqs= seq;% create a new array which is same as seq
            if (sum(seqs(1,roi)==1)/sum(seqs(1,roi)~=0))>thu %&& i~=144 && i~=145
                seqs(1,mot_up:mot_down)=1;
            end
            seqs(seqs<0) = 0;% replace all -1 in the array to 0
            score_all_pos(m,a) = sum(seqs(1,(a-border):(a+border)).*-score_matrix);%calculate the nucleosome dyad score for each position in the array. The score at each location is equal to the the sum of the sequence array mutiply by the score matrix
        end
        %% scan the read and place nucleosomes on the read
        for j = 1:bkg_nuc_number
            if j == 1 %set the default dyad of nuc-6, and the scanwindow
                if isempty(find(seq(1,1:100)==1))
                    pos_mov = (1+68+border):1:(lenmat-70-border);% If no methylated C in the first 100bp, set the border of the nuc at 69bp from the edge
                else
                    pos_mov = (1+border):1:(lenmat-70-border);%reset the moving window of the nucleosome
                end
            else       %set the default dyad of other nucs base on the position of the pre-defined upstream nucleosome
                if dyad_f(1,j-1)==0
                    break   %if no previous nucleosome, quit the cycle
                else
                    dyad_d = dyad_f(1,j-1)+2*b+15;% the initial pos calculated base on the downstream nuc
                    pos_mov = dyad_d:1:(lenmat-70-border);% set the scan window for this nucleosome
                    if lenmat-70-dyad_d<border
                        break % if no more space to put in an nucleosome, quit the cycle
                    end
                end
            end
            dyad=zeros(length(pos_mov),1);%create an array to save the dyad position of each window
            bound = zeros(length(pos_mov),2);%create a matrix to save the upstream and downstream boundary for the nucleosome
            for n = 1:length(pos_mov)
                dyad(n,1) = pos_mov(1,n);%set the dyad position in each window
                if (dyad(n,1)+border)>lenmat-70
                    break % if the edge of the nucleosome is out of bound, quit the cycle
                end
                bound(n,:) = [(dyad(n,1)-border),(dyad(n,1)+border)];% set the boundary based on the pos of dyad
                score_nuc = seqs(1,bound(n,1):bound(n,2)).*-score_matrix;% calculate the dyad score in each window
                %calculate the distance between adjacent
                %nucleosomes
                if j == 1
                    distance = 100;% for the first nucleosome ignore this distance
                else
                    distance = abs(dyad(n,1)-dyad_f(1,j-1))-2*b;% calculate the distance between this nuc and the upstream nuc
                end
                % calculate the penalty score based on the distance
                % between two nucleosomes
                pfactor = 1;%set the penalty factor
                penalty = pen_score(distance,pfactor,linker);% calculate penalty score using the function pen_score
                score_sum_new = sum(score_nuc) - penalty;% set the final score as the score - penalty score
                
                % if the nucleosome positioning score at the new
                % location is high than at the previous position, replace the old one with the new one.
                if score_sum_new >= -0.5 % set the minimium score for the existence of a nucleosome
                    if  score_sum_new > score_sum(j)
                        score_sum(j) = score_sum_new;%set the highest score
                        dyad_f(1,j) = dyad(n,1);%set the best dyad
                    else
                        break
                    end
                else
                end
            end%length(pos_mov)
        end
        %% calculate shifted nucleosome positions using peak_shift_COMB_motif_plc
        % the sequence, motif positions,predicted dyad positions, and nuc score were input
        [shift_all,shift_dyad,m1] = peak_shift_COMB_motif_plc(seqs,mpos, mot_up,mot_down,dyad_f,score_all_pos(m,:),b,c_limit,linker,bkg_nuc_number);
        nuc_shift_all(m,:) = shift_all;%
        shift_dyad_all(m,:)=shift_dyad;%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% adjust the position of nucleosome based on the dyads calculated above
        dyad_all(m,:) = dyad_f;
        pileup_m1(m,:) = m1;
        map(m,:) = seq;
    end%1-m
    % divid results into 2 structures, info4(pileup data) and infom(prediction and methylation matrices)
    pileup_met(1,:) = gc_pos;
    pileup_met(2,:) = mean(seqt(:,gc_pos));
    info4(i).pileup = mean(pileup_m1,1);
    info4(i).pileup_met = pileup_met;
    infom(i).pred = pileup_m1;
    infom(i).met = map;
end%1-169

%output data
save_data(pred_pileup,info4);
save_data(pred_matrix,infom);
disp('nucleosome predicted');
toc;
end

