function met_data_analysis_pip()
tic
%cd foldername;%open the folder that contains files downloaded from github
%%
mkdir('seqs');%make a folder named seqs.
ccs = 'ccs\example.fastq';%sequencing read file
seqs_rc = 'seqs\example.mat';%output file
seqs_c = rev2fw_header_v2(ccs,seqs_rc);%generate reverse complimentary reads

%%
mkdir('aligned');% make a directory named aligned
ref_bkg = 'ref\ref_bkg.txt';%the reference of the background sequence
ref_ndf = 'ref\ref_ndf.txt';%the reference of the sequence with an NDF motif
mot_bkg = 'ref\motif_bkg.mat';%the position of the variable region, the background sequence of the variable region
mot_ndf = 'ref\motif_ndf.mat';%the position of the motif, the sequence of the motif, and the name of the NDF
aligned_bkg_file = 'aligned\aligned_example_bkg.mat';%output file
aligned_ndf_file = 'aligned\aligned_example_ndf.mat';%output file
variable_region = 437:519;%the location of the variable region in the reference sequence


aligned_bkg = align_met_all_v2(seqs_c,ref_bkg,mot_bkg,variable_region,aligned_bkg_file);%align reads to the background sequence
aligned_ndf = align_met_all_v2(seqs_c,ref_ndf,mot_ndf,variable_region,aligned_ndf_file);%align reads to the NDF sequence

%%
mkdir('matrix');%make a directory named matrix
mat_bkg = 'matrix\matrix_example_bkg.mat';%output file for the bkg sequence
mat_ndf = 'matrix\matrix_example_ndf.mat';%output file for the NDF sequence
matrix_bkg = pbmatrix_v2(aligned_bkg,ref_bkg,mat_bkg);%generate a matrix for the bkg sequence
matrix_ndf = pbmatrix_v2(aligned_ndf,ref_ndf,mat_ndf);%generate a matrix for the NDF sequence


%%
border = 73;% the protection length on one side of the nucleosome. Full length is 73+1+73
bkg_pseudo_pos = [488, 508];%the position of pseudo motif on bkg sequence
is_remove_nuc = 1;%remove the nucleosome on ndf motifs if methylation level near motifs is high
is_bkg_yes = 1;%predict nucleosome on the bkg sequence
is_bkg_no = 0;%predict nucleosome on a ndf sequence
linker = 20;%No penalty for linkers longer than 20bp.
roi = 250:803;%region of interest(The HOpr region containing nucleosome -5, -4, and -3)
c_limit = 100;%The limit of optimization cycles
mkdir('nuc_prediction');%make directory called nuc_prediction
pred_matrix_bkg = 'nuc_prediction/nuc_pred_example_bkg.mat';% output bkg matrices
pred_pileup_bkg = 'nuc_prediction/nuc_pred_pileup_bkg.mat';% output averaged occupancies on the bkg sequence
pred_matrix_ndf = 'nuc_prediction/nuc_pred_example_ndf.mat';% output ndf matrices
pred_pileup_ndf = 'nuc_prediction/nuc_pred_pileup_ndf.mat';% output averaged occupancies on the ndf sequence


predicted_bkg = nuc_predict_v2(matrix_bkg,matrix_bkg,mot_ndf,bkg_pseudo_pos,border,linker,is_bkg_yes,is_remove_nuc,c_limit,roi,pred_matrix_bkg,pred_pileup_bkg);%predict nucleosomes on bkg reads
predicted_ndf = nuc_predict_v2(matrix_bkg,matrix_ndf,mot_ndf,bkg_pseudo_pos,border,linker,is_bkg_no,is_remove_nuc,c_limit,roi,pred_matrix_ndf,pred_pileup_ndf);%predict nucleosomes on ndf reads


%%

bkg_lp = 'nuc_prediction\pred_bkg_lp.mat';
ndf_lp = 'nuc_prediction\pred_ndf_lp.mat';
NDR_num_len_v2(predicted_bkg,mot_bkg,is_bkg_yes,bkg_pseudo_pos,border,bkg_lp);
NDR_num_len_v2(predicted_ndf,mot_ndf,is_bkg_no,bkg_pseudo_pos,border,ndf_lp);
%%
mkdir('nuc_prediction\heatmaps_pred');
mkdir('nuc_prediction\heatmaps_met');
hm_pred_bkg =  'nuc_prediction/heatmaps_pred/heatmap_pred_bkg';
hm_met_bkg = 'nuc_prediction/heatmaps_met/heatmap_met_bkg';
hm_pred_ndf =  'nuc_prediction/heatmaps_pred/heatmap_pred_ndf';
hm_met_ndf = 'nuc_prediction/heatmaps_met/heatmap_met_ndf';
is_sort_yes = 1;
simple_sort_yes = 1;
simple_sort_region = 150;
plot_region = 21:1174;
Pacbio_clustergram_v2(predicted_bkg,ref_bkg,mot_bkg,is_bkg_yes,bkg_pseudo_pos,plot_region,is_sort_yes,simple_sort_yes,simple_sort_region,hm_pred_bkg,hm_met_bkg);
Pacbio_clustergram_v2(predicted_ndf,ref_ndf,mot_ndf,is_bkg_no,bkg_pseudo_pos,plot_region,is_sort_yes,simple_sort_yes,simple_sort_region,hm_pred_ndf,hm_met_ndf);
toc

end