function Pacbio_clustergram_v2(predicted_bkg,predicted_tf,ref_bkg,ref_tf,mot_info,bkg_pseudo_pos,plot_region,is_sort,sort_region,heatmap_pred_bkg,heatmap_met_bkg,heatmap_pred_tf,heatmap_met_tf)

Pacbio_clustergram(predicted_bkg,ref_bkg,mot_info,1,bkg_pseudo_pos,plot_region,is_sort,sort_region,heatmap_pred_bkg,heatmap_met_bkg);%make heatmaps for the background sequence
Pacbio_clustergram(predicted_tf,ref_tf,mot_info,0,bkg_pseudo_pos,plot_region,is_sort,sort_region,heatmap_pred_tf,heatmap_met_tf);%make heatmaps for sequences with TF motifs
end

function Pacbio_clustergram(predicted,refs,mot,is_bkg,bkg_pseudo_pos,plot_region,is_sort,sort_region,hm_pred,hm_met)
%% This function plot heatmaps using methylation and nucleosome prediction matrices
%predicted is the output of nuc_predict_v2.m. It contains the nucleosome prediction matrix.
%refs is the reference file which is a txt file.
%mot is the motif information file which is a mat file containing the position, sequence, and name of each TF motif
%is_bkg determines whether the seqence is a background sequence. 1 is yes. 0 is no.
%bkg_pseudo_pos is the pseudo position of the background sequence
%plot_region is the region in the matrix that will be plotted.
%is_sort determines whether rows in the heatmap are sorted. 1 is yes. 0 is no.
%sort_region is the region near motifs that will be used to sort heatmaps
%hm_pred is the output for nucleosome prediction map
%hm_met is the output for methylation protection map

references = readcell(refs);
%% plot heatmaps
dbstop if error
%% plot heatmaps for each read
for m = 1:length(references)
    seqs = predicted(m).met; % read the methylation info
    nuc = predicted(m).pred; % read the nucleosome info
    ref = references{m,1};

    if is_bkg ==1%detemine whether the sequence has motifs or not.
        mpos = bkg_pseudo_pos;
        mname = 'bkg';
    else
        mpos = load(mot).mot.tf(m).pos;%read motif positions
        mname = load(mot).mot.tf(m).name;
    end

    gc_pos = strfind(ref,'GY');
    map = protection_strict_map(seqs,gc_pos(1,1),gc_pos(end,end)); % make a protection map based on the position of GCs/GTs
    [len,~] = size(nuc);
    if is_sort == 1
        sp_score = sum(nuc(:,(mpos(1,1)-sort_region):(mpos(end,end)+sort_region)),2);
        nuc_sp = [transpose(1:len) sp_score];
        [~,idx] = sort(nuc_sp(:,2),'ascend');
        sort_sp = nuc_sp(idx,:);
        sort_all = [sort_sp nuc(idx,:)];
        sort_map = map(idx,:);
    else
        sort_all = nuc;
        sort_map = seqs;
    end
    %% change all values to 1 or 0 for better visualization, point out the motifs
    sort_all(sort_all>0)=1;
    sort_all(sort_all<=0)=0;
    sort_all(:,[mpos(1,1),mpos(end,end)]-2)= 2;
    sort_map(:,[mpos(1,1),mpos(end,end)])=2;

    %% plot and save heatmaps
    fig1 = heatmap(1:length(plot_region),1:len,sort_all(:,(plot_region+2)),'GridVisible','off','ColorbarVisible','off','Colormap',[1 1 1; 0.29 0.66 0.78;0 0 0]);
    print([hm_pred,'_',mname,'_',num2str(m),'.png'],'-dpng','-r300');
    savefig([hm_pred,'_',mname,'_',num2str(m),'.fig']);
    close all force
    fig2 = heatmap(1:length(plot_region),1:len,sort_map(:,plot_region),'GridVisible','off','ColorbarVisible','off','Colormap',[0.996 0.6 0.6; 1 1 1; 0.29 0.66 0.78; 0 0 0]);
    print([hm_met,'_',mname,'_',num2str(m),'.png'],'-dpng','-r300');
    savefig([hm_met,'_',mname,'_',num2str(m),'.fig']);
    close all force
end
end
