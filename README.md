# Mapping nucleosome occupancy 
We present a method to map nucleosomes using a DNA methyltransferase. The scripts here are used to align bisulfite-converted long Pacbio reads and predict nucleosomes on these reads. For complete details on the use and execution of this method, please refer to "[Chen et al. 2022](https://www.cell.com/cell-reports/fulltext/S2211-1247(22)01068-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124722010683%3Fshowall%3Dtrue)"
Bioinformatics and curve fitting tool boxes are required for data analysis.

## Input files
### Pacbio sequencing reads
Reads are in fastq format. The example file "[example.fastq](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/consensus_sequence/example.fastq)" contains 6000 reads. The method for generating these reads can be found in "[Chen et al. 2022](https://www.cell.com/cell-reports/fulltext/S2211-1247(22)01068-3?_returnURL=https%3A%2F%2Flinkinghub.elsevier.com%2Fretrieve%2Fpii%2FS2211124722010683%3Fshowall%3Dtrue)".

### Reference sequences
Reference sequences are in txt format. The example file "[ref_bkg.txt](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/ref/ref_bkg.txt)" and "[ref_tf.txt](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/ref/ref_tf.txt)" contains the background sequence and sequences with TF motifs respectively.
### Positions and names of motifs
The motif information is saved in mat format. The example file "[motif_info.mat](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/ref/motif_info.mat)" contains the motif information of TF references and the sequence information of the background sequence.
### MATLAB scripts
All MATLAB scripts can be found in the main folder(https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo).

## Nucleosome prediction
In MATLAB script, words following % are comments.
### Generate reverse complimentary reads
Run the command below to generate reverse complimentary reads of the raw reads
> cd yourFolder;%open the folder containing all files. Replace "yourFolder" with your directory path.\
> mkdir('seqs');%make a folder named seqs.\
> ccs_file = 'consensus_sequence\example.fastq';%sequencing read file\
> seqs_file = 'seqs\example.mat';%output file which contains both forward and reverse complimentary reads.\
> seqs_all = rev2fw_header_v2(ccs_file,seqs_file);%generate reverse complimentary reads

Forward and reverse reads are saved in a MATLAB structure. 

### Align reads to reference sequences
Run the command below to align reads to reference sequences
> mkdir('aligned');% make a directory named aligned\
> ref_bkg = 'ref\ref_bkg.txt';%the reference of the background sequence\
> ref_tf = 'ref\ref_tf.txt';%the reference of the sequence with an TF motif\
> mot_info = 'ref\motif_info.mat';%the position of the motif, the sequence of the motif, and the name of the TF on each reference sequence\
> aligned_bkg_file = 'aligned\aligned_example_bkg.mat';%output file\
> aligned_tf_file = 'aligned\aligned_example_tf.mat';%output file\
> variable_region = 437:519;%the location of the variable region in the reference sequence\
> align_accuracy = 0.98;%The minimal alignment accuracy is 98%\
> is_bkg = 1;% align reads to the bkg sequence\
> aligned_bkg = align_met_all_v2(seqs_all,is_bkg,ref_bkg,mot_info,variable_region,align_accuracy,aligned_bkg_file);%align reads to the background sequence\
> is_bkg = 0;% align reads to sequences containing TF motifs\
> aligned_tf = align_met_all_v2(seqs_all,is_bkg,ref_tf,mot_info,variable_region,align_accuracy,aligned_tf_file);%align reads to the TF sequence

align_met_all_v2.m aligned long reads to reference sequences using global and local alignment functions in the MATLAB bioinformatics toolbox.
In this example, after this alignment step, all aligned reads are adjusted to 1278bp.

### Generate methylation matrixes
Run the command below to convert aligned reads to numeric arrays. Methylated GCs are 1, while other nucleotides are 0.
> mkdir('matrix');%make a directory named matrix\
> matrix_bkg_file = 'matrix\matrix_example_bkg.mat';%output file for the bkg sequence\
> matrix_tf_file = 'matrix\matrix_example_tf.mat';%output file for the TF sequence\
> matrix_bkg = pbmatrix_v2(aligned_bkg,ref_bkg,matrix_bkg_file);%generate a matrix for the bkg sequence\
> matrix_tf = pbmatrix_v2(aligned_tf,ref_tf,matrix_tf_file);%generate a matrix for the TF sequence

pbmatrix_v2.m convert DNA sequences to numeric arrays. GCs are converted to 1 while other nucleotides are converted to 0.
All reads aligned to one reference sequence will be combined into a N by 1278 matrix. N is the number of reads.
The output file is in MATLAB structure format. Matrixes are saved separately in the structure.

### Predict nucleosomes
Run the command below to predict nucleosomes on each sequencing reads.
> border = 73;% the protection length on one side of the nucleosome. Full length is 73+1+73\
> bkg_pseudo_pos = [488, 508];%the position of pseudo motif on bkg sequence\
> is_remove_nuc = 1;%remove the nucleosome on TF motifs if methylation level near motifs is high\
> linker = 20;%No penalty for linkers longer than 20bp.\
> roi = 250:803;%region of interest(The HOpr region containing nucleosome -5, -4, and -3)\
> c_limit = 100;%The limit of optimization cycles\
> mkdir('nuc_prediction');%make directory called nuc_prediction\
> pred_matrix_bkg = 'nuc_prediction/nuc_pred_example_bkg.mat';% output bkg matrices\
> pred_pileup_bkg = 'nuc_prediction/nuc_pred_pileup_bkg.mat';% output averaged occupancies on the bkg sequence\
> pred_matrix_tf = 'nuc_prediction/nuc_pred_example_tf.mat';% output TF matrices\
> pred_pileup_tf = 'nuc_prediction/nuc_pred_pileup_tf.mat';% output averaged occupancies on the TF sequence\
> is_bkg = 1;%predict nucleosome on the bkg sequence\
> predicted_bkg = nuc_predict_v2(matrix_bkg,matrix_bkg,mot_info,bkg_pseudo_pos,border,linker,is_bkg,is_remove_nuc,c_limit,roi,pred_matrix_bkg,pred_pileup_bkg);%predict nucleosomes on bkg reads\
> is_bkg = 0;%predict nucleosome on a TF sequence\
> predicted_tf = nuc_predict_v2(matrix_bkg,matrix_tf,mot_info,bkg_pseudo_pos,border,linker,is_bkg,is_remove_nuc,c_limit,roi,pred_matrix_tf,pred_pileup_tf);%predict nucleosomes on TF reads

nuc_predict_v2.m predicts nucleosomes using matrixes generated by pbmatrix_v2.m. It generates a N by 1278 matrix for each reference sequence. In the matrix, one row represent the nucleosome occupancy on a sequence.
The output file is in MATLAB structure format.

### Calculate NDR length and probability
Run the command below to calculate NDR length and probability from the nucleosome prediction matrixes
> bkg_lp = 'nuc_prediction\pred_bkg_lp.mat';%output file\
> tf_lp = 'nuc_prediction\pred_tf_lp.mat';%output file\
> is_bkg = 1;%predict nucleosome on the bkg sequence\
> NDR_num_len_v2(predicted_bkg,mot_info,is_bkg,bkg_pseudo_pos,border,bkg_lp);%calculate length and probability of NDRs on bkg reads\
> is_bkg = 0;%predict nucleosome on a TF sequence\
> NDR_num_len_v2(predicted_tf,mot_info,is_bkg,bkg_pseudo_pos,border,tf_lp);%calculate length and probability of NDRs on TF reads

NDR_num_len_v2.m calculate the length and proportion of NDRs. 
The output file is in MATLAB structure format.

### Make methylation and nucleosome prediction heatmaps
Run the command below to generate methylation and nucleosome prediction heatmaps. 
> mkdir('nuc_prediction\heatmaps_pred');%make a directory named heatmaps_pred\
> mkdir('nuc_prediction\heatmaps_met');%make a directory named heatmaps_met\
> heatmap_pred_bkg =  'nuc_prediction/heatmaps_pred/heatmap_pred_bkg';%output nucleosome prediction heatmap for the background sequence\
> heatmap_met_bkg = 'nuc_prediction/heatmaps_met/heatmap_met_bkg';%output methylation map for the background sequence\
> heatmap_pred_tf =  'nuc_prediction/heatmaps_pred/heatmap_pred_tf';%output nucleosome prediction maps for the TF sequences\
> heatmap_met_tf = 'nuc_prediction/heatmaps_met/heatmap_met_tf';%output methylation maps for the TF sequences\
> is_sort = 1;%sort reads when make heatmaps\
> sort_region = 150;%sort reads based on the ±150bp near motifs\
> plot_region = 21:1174;%plot 21 to 1174 columns in the matrix\
> is_bkg = 1;%make the heatmap for the bkg sequence\
> Pacbio_clustergram_v2(predicted_bkg,ref_bkg,mot_info,is_bkg,bkg_pseudo_pos,plot_region,is_sort,sort_region,heatmap_pred_bkg,heatmap_met_bkg);%make heatmaps\
> is_bkg = 0;%make heatmaps for TF sequences\
> Pacbio_clustergram_v2(predicted_tf,ref_tf,mot_info,is_bkg,bkg_pseudo_pos,plot_region,is_sort,sort_region,heatmap_pred_tf,heatmap_met_tf);%make heatmaps

Pacbio_clustergram_v2.m generate heatmaps using the matrixes made by pbmatrix_v2.m and nuc_predict_v2.m. Plots are saved in MATLAB fig and jpg formats.

## Output files
All example output files can be found in folder "[example_output](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/tree/main/example_output)". You can compare the output files you generated with the example files to make sure the scripts work properly.
### forward and reverse complimentary reads
Forward and reverse complimentary reads are saved in "[example.mat](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/example_output/seqs/example.mat)".

### aligned reads
Aligned reads are saved in folder "[aligned](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/tree/main/example_output/aligned)".

### Methylation matrix
Methylation matrixes are saved in folder "[matrix](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/tree/main/example_output/matrix)".

### Nucleosome occupancy matrix
Nucleosome occupancy matrixes are saved in "[nuc_pred_example_bkg.mat](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/example_output/nuc_prediction/nuc_pred_example_bkg.mat)" and "[nuc_pred_example_tf.mat](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/example_output/nuc_prediction/nuc_pred_example_tf.mat)".

### NDR length and probability
NDR length and probability are saved in "[pred_bkg_lp.mat.mat](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/example_output/nuc_prediction/pred_bkg_lp.mat.mat)" and "[pred_tf_lp.mat.mat](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/blob/main/example_output/nuc_prediction/pred_tf_lp.mat.mat)".

### Methylation and nucleosome prediction bheatmaps
Output heatmaps are saved in folders "[heatmaps_met](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/tree/main/example_output/nuc_prediction/heatmaps_met)" and "[heatmaps_pred](https://github.com/HengyeChen/Pacbio_nucleosome_mapping_demo/tree/main/example_output/nuc_prediction/heatmaps_pred)". 
