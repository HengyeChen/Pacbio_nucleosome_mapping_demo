function met_data_analysis_pip(name,r_index)
tic
name = 'example';
r_index = [1, 2, 144, 145];%144 and 145 are background sequences

mkdir('ccs');
mkdir('ref');


rev2fw_header_v2(name);
disp('rev2fw_header_v2 completed');

align_met_all_v2(name,r_index);
disp('alignment completed');

pbmatrix_v2(name,r_index);
disp('matrix generated');

nuc_predict_v2(name,r_index);
disp('nucleosome predicted');

NDR_num_len_v2(name,r_index);
disp('NDR length and proportion calculated');

Pacbio_clustergram_v2(name,r_index);
disp('heatmaps plotted');
toc

end