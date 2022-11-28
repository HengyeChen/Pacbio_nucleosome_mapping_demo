function [info4, infom] = nuc_predict_v2(name,r_index)%This version predict nuc-6 at the canonical place, when no methylated GC at the upstream of the sequence
%% Parameter initiation

mkdir('nuc_prediction');
motifs = load('ref/motif_pos_v2.mat');%read motif positions
mkdir('nuc_prediction/tt');%make folder to save nucleosome prediction data

%%%%%%%%%%%%%%%%%% generate a cumulative distribution curve, and use it as the scoring matrix for a mononucleosome
border = 73;% set the distance between the boundary and the dyad
linker = 20;% set the length of standard linker. When predicted linker is short than this value, set a penalty for it. The minmal length of linker is 5bp shorter than this value.
a = -20; b = border;
x = a:1:b;%set the range of the cumulative distribution curve
m = (a + b)/2;%calculate the mean of a and b, which also is the middle point of x
s = 10; % set the standard deviation of the curve
pd = makedist('Normal','mu',m,'sigma',s);%creates a probability Normal distribution object in which 'mu'(the mean) is m, and 'sigma'(standard deviation) is s
f = cdf(pd,x);%Create a standard normal distribution object with the mean, ? is m and the standard deviation, ? is s. Each value in f correspond to a value in x.
score_matrix = [f(1,22:end), 1, flip(f(1,22:end))];%trim the distribution array, flip the trimmed array, and combine these two arrays. The final array is the score matrix for a nulceosome

c_limit = 100; % cycle limit for nucleosome positioning refinement
roi = 250:803;
%% nucleosome prediction
    %name = list{k,1};% select the sample name from the list
    data = load(['matrix\matrix_',name,'.mat']);
    tic;
    dbstop if error
    thub= [0,0];
    gcn = [19,17];
    %% set the methylation threshold with the methylation level in background sequences
    for i = 144:145 %144 and 145 are background sequence
        met_all = data.C_T_all_sum(i).C_T_sum_trim;
        mm = met_all(:,roi);
        mm(mm==-1)=0;
        xxx = sum(mm,2)/gcn(1,i-143);
        ds = fitdist(xxx,'Normal');
        %thub(1,i-143) = ds.mu + 1.645*ds.sigma;
        thub(1,i-143) = ds.mu + 1.282*ds.sigma;
    end
    thu = mean(thub);
    %% predict nucleosome in each sequence
    for ri = 1:length(r_index)
        i = r_index(ri);
        mpos = motifs.data(i).pos;%read the position file
        if i == 144 || i == 145
            mpos = [495;505];
        end
        matrix = data.C_T_all_sum(i).C_T_sum_trim;%read the trimmed sequence matrix
        [len,~] = size(matrix);%count the number of reads in the matrix
        if len == 0
            continue
        end
        dyad_all = zeros(len,7);%create an zero matrix for the positions of dyads of all sequences
        %% predict nucleosome positions
        map = zeros(len,1278);%
        nuc_shift_all = zeros(len,1278);% create an zero matrix for the shifted nucleosome positions
        score_all_pos = zeros(len,1278);
        shift_dyad_all = zeros(len,7);% create an zero matrix for the shifted dyad positions
        mot_down = max(max(mpos))+10;%set the boundary of the motifs, downstream + 10bp
        mot_up = min(min(mpos))-10;%set the boundary of the motifs, upstream - 10bp
        pileup_m1 = zeros(len,1278);
        gc_pos = readmatrix('ref\GC_positions.xlsx','Sheet',i);
        pileup_met = zeros(3,length(gc_pos));
        seqt = matrix;% create a new array which is same as seq
        seqt(seqt<0) = 0;
        for m = 1:len
            seq = matrix(m,:);%read a single sequence from the matrix
            dyad_f = zeros(1,7);%reset dyad positions of all nucleosomes
            score_sum = zeros(1,7)-2;% set the default score of each nucleosome as -2
            
            for a = (1+border):(1278-border)
                seqs= seq;% create a new array which is same as seq
                if (sum(seqs(1,roi)==1)/sum(seqs(1,roi)~=0))>thu %&& i~=144 && i~=145
                    seqs(1,mot_up:mot_down)=1;
                end
                seqs(seqs<0) = 0;% replace all -1 in the array to 0
                score_all_pos(m,a) = sum(seqs(1,(a-border):(a+border)).*-score_matrix);%calculate the nucleosome dyad score for each position in the array. The score at each location is equal to the the sum of the sequence array mutiply by the score matrix
            end
            %% scan the read and place nucleosomes on the read
            for j = 1:7
                if j == 1 %set the default dyad of nuc-6, and the scanwindow
                    if isempty(find(seq(1,1:100)==1))
                        pos_mov = (1+68+border):1:(1208-border);% If no methylated C in the first 100bp, set the border of the nuc at 69bp from the edge
                    else
                        pos_mov = (1+border):1:(1208-border);%reset the moving window of the nucleosome
                    end
                else       %set the default dyad of other nucs base on the position of the pre-defined upstream nucleosome
                    if dyad_f(1,j-1)==0
                        break   %if no previous nucleosome, quit the cycle
                    else
                        dyad_d = dyad_f(1,j-1)+2*b+15;% the initial pos calculated base on the downstream nuc
                        pos_mov = dyad_d:1:(1208-border);% set the scan window for this nucleosome
                        if 1208-dyad_d<border
                            break % if no more space to put in an nucleosome, quit the cycle
                        end
                    end
                end
                dyad=zeros(length(pos_mov),1);%create an array to save the dyad position of each window
                bound = zeros(length(pos_mov),2);%create a matrix to save the upstream and downstream boundary for the nucleosome
                for n = 1:length(pos_mov)
                    dyad(n,1) = pos_mov(1,n);%set the dyad position in each window
                    if (dyad(n,1)+border)>1208
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
            [shift_all,shift_dyad,m1] = peak_shift_COMB_motif_plc(seqs,mpos, mot_up,mot_down,dyad_f,score_all_pos(m,:),i,b,c_limit,linker);
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
    toc;
    %output data
   save_data(['nuc_prediction/tt/nuc_pred_pileup_all_TTxx_',name,'.mat'],info4);
   save_data(['nuc_prediction/tt/nuc_pred_TTxx_',name,'.mat'],infom);
   
end

