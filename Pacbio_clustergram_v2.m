function Pacbio_clustergram_v2(predicted,refs,mot,is_bkg,bkg_pseudo_pos,plot_region,is_sort,simple_sort,simple_sort_region,hm_pred,hm_met)
%% This function plot heatmaps using methylation and nucleosome prediction matrices
%predicted is the output of nuc_predict_v2.m. It contains the nucleosome prediction matrix.
%refs is the reference file which is a txt file.
%mot is the motif information file which is a mat file containing the position, sequence, and name of each TF motif
%is_bkg determines whether the seqence is a background sequence. 1 is yes. 0 is no.
%bkg_pseudo_pos is the pseudo position of the background sequence
%plot_region is the region in the matrix that will be plotted.
%is_sort determines whether rows in the heatmap are sorted. 1 is yes. 0 is no.
%simple_sort determines whether rows in the heatmap are sorted in the simple way based on the postion of motifs. 1 is yes. 0 is no.
%simple_sort_region is the region near motifs that will be used to sort heatmaps
%hm_pred is the output for nucleosome prediction map
%hm_met is the output for methylation protection map 



%%detemine whether the sequence has motifs or not.
if is_bkg ==1
    mpos = bkg_pseudo_pos;
else
    mpos = load(mot).mot.pos;%read motif positions
end
references = readcell(refs);
%% plot heatmaps
dbstop if error
%% plot heatmaps for each read
for m = 1
    seqs = predicted.met; % read the methylation info
    nuc = predicted.pred; % read the nucleosome info
    ref = references{m,1};
    gc_pos = strfind(ref,'GY');
    map = protection_strict_map(seqs,gc_pos(1,1),gc_pos(end,end)); % make a protection map based on the position of GCs/GTs
    [len,~] = size(nuc);
    if is_sort == 1 && simple_sort == 1
        sp_score = sum(nuc(:,(mpos(1,1)-simple_sort_region):(mpos(end,end)+simple_sort_region)),2);
        nuc_sp = [transpose(1:len) sp_score];
        [~,idx] = sort(nuc_sp(:,2),'ascend');
        sort_sp = nuc_sp(idx,:);
        sort_all = [sort_sp nuc(idx,:)];
        sort_map = map(idx,:);
    elseif is_sort == 1 && simple_sort == 0
        %% calculate the methylation in different regions
        score4 = sum(nuc(:,400:550),2);
        score3 = sum(nuc(:,551:750),2);
        score5 = sum(nuc(:,250:399),2);
        score21 = sum(nuc(:,751:1050),2);
        score6 = sum(nuc(:,100:249),2);
        %% calculate a methylation score for each read based on the methylation at different regions and sort reads by this score
        nuc_43 = [transpose(1:len) (score5*5 + score4*100 + score3*10+score21*2+score6*0.5)];
        [~,idx] = sort(nuc_43(:,2),'ascend');
        sort_43 = nuc_43(idx,:);
        sort_all = [sort_43 nuc(idx,:)];
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
end

%% plot and save heatmaps
fig1 = heatmap(1:length(plot_region),1:len,sort_all(:,(plot_region+2)),'GridVisible','off','ColorbarVisible','off','Colormap',[1 1 1; 0.29 0.66 0.78;0 0 0]);
print([hm_pred,'.png'],'-dpng','-r300');
savefig([hm_pred,'.fig']);
close all force
fig2 = heatmap(1:length(plot_region),1:len,sort_map(:,plot_region),'GridVisible','off','ColorbarVisible','off','Colormap',[0.996 0.6 0.6; 1 1 1; 0.29 0.66 0.78; 0 0 0]);
print([hm_met,'.png'],'-dpng','-r300');
savefig([hm_met,'.fig']);
close all force
end
