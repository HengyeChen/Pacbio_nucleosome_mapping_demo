function met_data_analysis_pip()
tic
%open the folder that contains files downloaded from github
%cd yourFolder;

%% generate reverse complimentary reads of consensus reads
mkdir('seqs');%make a folder named seqs.
ccs_file = 'consensus_sequence\example.fastq';%sequencing read file
seqs_file = 'seqs\seqs_all_record.mat';%output file
seqs_all = rev2fw_header_v2(ccs_file,seqs_file);%generate reverse complimentary reads

%% align reads to references
mkdir('aligned');% make a directory named aligned
ref_bkg = 'ref\ref_bkg.txt';%the reference of the background sequence
ref_tf = 'ref\ref_tf.txt';%the reference of the sequence with an TF motif
mot_info = 'ref\motif_info.mat';%the position of the motif, the sequence of the motif, and the name of the TF on each reference sequence
aligned_bkg_file = 'aligned\aligned_example_bkg.mat';%output file
aligned_tf_file = 'aligned\aligned_example_tf.mat';%output file
variable_region = 437:519;%the location of the variable region in the reference sequence
align_accuracy = 0.98;%The minimal alignment accuracy is 98%

[aligned_bkg, aligned_tf] = align_met_all_v2(seqs_all,ref_bkg,ref_tf,mot_info,variable_region,align_accuracy,aligned_bkg_file,aligned_tf_file);%align reads to the reference sequences

%% convert aligned reads to numeric matrices
mkdir('matrix');%make a directory named matrix
matrix_bkg_file = 'matrix\matrix_example_bkg.mat';%output file for the bkg sequence
matrix_tf_file = 'matrix\matrix_example_tf.mat';%output file for the TF sequence
matrix_bkg = pbmatrix_v2(aligned_bkg,ref_bkg,matrix_bkg_file);%generate a matrix for the bkg sequence
matrix_tf = pbmatrix_v2(aligned_tf,ref_tf,matrix_tf_file);%generate a matrix for the TF sequence

%% predict nucleosomes
border = 73;% the protection length on one side of the nucleosome. Full length is 73+1+73
bkg_pseudo_pos = [488, 508];%the position of pseudo motif on bkg sequence
is_remove_nuc = 1;%remove the nucleosome on TF motifs if methylation level near motifs is high
linker = 20;%No penalty for linkers longer than 20bp.
roi = 250:803;%region of interest(The HOpr region containing nucleosome -5, -4, and -3)
c_limit = 100;%The limit of optimization cycles
mkdir('nuc_prediction');%make directory called nuc_prediction
pred_matrix_bkg = 'nuc_prediction/nuc_pred_example_bkg.mat';% output bkg matrices
pred_pileup_bkg = 'nuc_prediction/nuc_pred_pileup_bkg.mat';% output averaged occupancies on the bkg sequence
pred_matrix_tf = 'nuc_prediction/nuc_pred_example_tf.mat';% output TF matrices
pred_pileup_tf = 'nuc_prediction/nuc_pred_pileup_tf.mat';% output averaged occupancies on the TF sequence

[predicted_bkg,predicted_tf] = nuc_predict_v2(matrix_bkg,matrix_tf,mot_info,bkg_pseudo_pos,border,linker,is_remove_nuc,c_limit,roi,pred_matrix_bkg,pred_pileup_bkg,pred_matrix_tf,pred_pileup_tf);%predict nucleosomes on reads

%% calculate NDR length and proportion
bkg_lp = 'nuc_prediction\pred_bkg_lp.mat';%output file
tf_lp = 'nuc_prediction\pred_tf_lp.mat';%output file

NDR_num_len_v2(predicted_bkg,predicted_tf,mot_info,bkg_pseudo_pos,border,bkg_lp,tf_lp);%calculate length and probability of NDRs on reads

%% make heatmaps
mkdir('nuc_prediction\heatmaps_pred');%make a directory named heatmaps_pred
mkdir('nuc_prediction\heatmaps_met');%make a directory named heatmaps_met
heatmap_pred_bkg =  'nuc_prediction/heatmaps_pred/heatmap_pred_bkg';%output nucleosome prediction heatmap for the background sequence
heatmap_met_bkg = 'nuc_prediction/heatmaps_met/heatmap_met_bkg';%output methylation map for the background sequence
heatmap_pred_tf =  'nuc_prediction/heatmaps_pred/heatmap_pred_tf';%output nucleosome prediction maps for the TF sequences
heatmap_met_tf = 'nuc_prediction/heatmaps_met/heatmap_met_tf';%output methylation maps for the TF sequences
is_sort = 1;%sort reads when make heatmaps
sort_region = 150;%sort reads based on the ±150bp near motifs
plot_region = 21:1174;%plot 21 to 1174 columns in the matrix

Pacbio_clustergram_v2(predicted_bkg,predicted_tf,ref_bkg,ref_tf,mot_info,bkg_pseudo_pos,plot_region,is_sort,sort_region,heatmap_pred_bkg,heatmap_met_bkg,heatmap_pred_tf,heatmap_met_tf);%make heatmaps
toc
end
