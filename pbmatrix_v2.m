function  pbmatrix_v2(name,r_index)
%% this function convert DNA sequences to numeric arrays

%list = [{'17R'};{'17D'};{'ID'}; {'IR'};{'ND'};{'NR'};{'TD'};{'TR'}];% set the name of samples

mkdir('matrix');

datam = load(['aligned\aligned_',name,'.mat']);

seqs_m_all = datam.alignment;
C_T_all_sum = struct('C_T_sum_trim',[],'C_T_pos',[],'value_sum_trim',[]);
for ki = 1:length(r_index)
    k = r_index(ki);
    seqs_c = seqs_m_all(1,k).reads;
    
    position = xlsread('ref\GC_positions.xlsx', k);
    [~, num_position] = size(position);
    value_sum = zeros(3,num_position);
    value_sum(1,:)=position;
    value_sum(3,:)=position-1208;
    [lenc, ~] = size(seqs_c);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [C_T_sum_c,~]= countreads(position,seqs_c,value_sum);
    C_T_sum_trim = zeros(lenc,1278);
    j=0;
    for i = 1:lenc
        if sum(C_T_sum_c(i,:)==1)< 4
            
        else
            j=j+1;
            C_T_sum_trim(j,:) = C_T_sum_c(i,:);
        end
    end
    C_T_sum_trim((j+1):end,:) = [];
    
    [~,widt] = size(C_T_sum_trim);
    value_sum_trim = sum(C_T_sum_trim)/widt;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    C_T_all_sum(k).C_T_sum_trim = C_T_sum_trim;
    C_T_all_sum(k).C_T_pos = position(1,:);%write C_T_sum and position to a structure
    C_T_all_sum(k).value_sum_trim = value_sum_trim;

end
save(['matrix\matrix_',name,'.mat'],'C_T_all_sum');

i =1;
