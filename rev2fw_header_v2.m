function seqs_c = rev2fw_header_v2(ccs_fq,seqs_rc)
%%
%seqs_c contains raw reads and reverse complementary reads
%ccs_fq is the input file which should be a fastq file
%seqs_rc is the output file

tic
[header,seqs] = fastqread(ccs_fq);% read the headers and sequences from the fastq file
%% generate a new file containing reads and reversed reads
[~,len] = size(seqs);
seqs_c = cell(len*2,2);% create a new cell from headers and sequences
for i = 1:len%input headers, seqs, and the reverse complementary seqs to seqs_c
    seq = seqs{1,i};% get one read from the seqs
    seqs_c{i,2} = header{1,i};% set the second column of seqs_c as headers
    seqs_c{len+i,2} = header{1,i};
    seqs_c{i,1}= seq;
    seqs_c{len+i,1}= seqrcomplement(seq);
end
save_data(seqs_rc,seqs_c);
disp('rev2fw_header_v2 completed');
toc
end
