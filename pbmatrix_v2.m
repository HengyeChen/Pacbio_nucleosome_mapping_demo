function  C_T_all_sum = pbmatrix_v2(aligned,refs,mat_output)
%% this function convert DNA sequences to numeric arrays
%aligned is the output of align_met_all_v2.m. It contains aligned reads 
%refs is the reference file which is a txt file.
%mat_output is the output file

tic
references = readcell(refs);%open the reference file for bkg sequences(orignial Cs are converted to Ys)
C_T_all_sum = struct('C_T_sum_trim',[],'C_T_pos',[],'value_sum_trim',[]);
%%
for k = 1
    seqs = aligned.reads;
    ref = references{k};
    position = strfind(ref,'GY');
    [~, num_position] = size(position);
    value_sum = zeros(3,num_position);
    value_sum(1,:)=position;
    value_sum(3,:)=position-1208;
    [lenc, ~] = size(seqs);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [C_T_sum_c,~]= countreads(position,seqs,value_sum);
    C_T_sum_trim = zeros(lenc,1278);
    j=0;
    %%reads with less than 4 GCs are removed
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
save(mat_output,'C_T_all_sum');
disp('matrix generated');
toc
end
