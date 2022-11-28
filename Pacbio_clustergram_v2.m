function Pacbio_clustergram_v2(name,r_index)
%% initiate parameters
%fname = '/storage/home/h/hxc62/work/';%define the home folder
%addpath(fname);
%fname = 'E:/d/Dropbox/';%set the home folder
%fname = 'C:\Users\bailab\Dropbox\';
%list = [{'NR3C'};{'OR3C'};{'TR'};{'WR2C'};{'IR2C'};{'WTcomb'}];% set the names of samples
%list = [{'IR2C'};{'NR3C'};{'OR3C'};{'TR'};{'WTall'};{'OR3Cn'};{'NR3Cn'};{'TRall'}];% set the names of samples
motifs = load('ref/motif_pos_v2.mat');%read motif positions
mkdir('nuc_prediction\heatmaps_pred');
mkdir('nuc_prediction\heatmaps_met');
%% plot heatmaps
%% read nucleosome and methylation info
%name = list{n,1};% select the sample name from the list
data = load(['nuc_prediction\tt\nuc_pred_TTxx_',name,'.mat']);
mkdir(['nuc_prediction\heatmaps_pred\pred_HC_cluster_colorc_TTxx_',name,'_mt']);
mkdir(['nuc_prediction\heatmaps_met\pred_HC_cluster_colorc_TTxx_',name,'_mt']);
dbstop if error
mm = r_index;% select sequences of interest
%% plot heatmaps for each read
for x = 1:length(mm)
    m = mm(1,x);% get the index
    mpos = motifs.data(m).pos;%read the position file
    seqs = data.data(m).met; % read the methylation info
    nuc = data.data(m).pred; % read the nucleosome info
    map = protection_strict_map(seqs); % make a protection map based on the position of GCs/GTs
    [len,~] = size(nuc);
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
    %% change all values to 1 or 0 for better visualization, point out the motifs
    sort_all(sort_all>0)=1;
    sort_all(sort_all<=0)=0;
    sort_all(:,[mpos(1,1),mpos(end,end)]-2)= 2;
    sort_map(:,[mpos(1,1),mpos(end,end)])=2;
    %% plot and save heatmaps
    fig1 = heatmap(1:1154,1:len,sort_all(:,23:1176),'GridVisible','off','ColorbarVisible','off','Colormap',[1 1 1; 0.29 0.66 0.78;0 0 0]);
    print(['nuc_prediction\heatmaps_pred\pred_HC_cluster_colorc_TTxx_',name,'_mt\pred_HC_cluster_mndr_',name,'_',num2str(m),'.png'],'-dpng','-r300');
    savefig(['nuc_prediction\heatmaps_pred\pred_HC_cluster_colorc_TTxx_',name,'_mt\pred_HC_cluster_mndr_',name,'_',num2str(m),'.fig']);
    close all force
    fig2 = heatmap(1:1154,1:len,sort_map(:,21:1174),'GridVisible','off','ColorbarVisible','off','Colormap',[0.996 0.6 0.6; 1 1 1; 0.29 0.66 0.78; 0 0 0]);
    print(['nuc_prediction\heatmaps_met\pred_HC_cluster_colorc_TTxx_',name,'_mt\pred_HC_cluster_mndr_',name,'_',num2str(m),'_GCmap.png'],'-dpng','-r300');
    savefig(['nuc_prediction\heatmaps_met\pred_HC_cluster_colorc_TTxx_',name,'_mt\pred_HC_cluster_mndr_',name,'_',num2str(m),'_GCmap.fig']);
    close all force
end
end