function rev2fw_header_v2(name)

%list = [{'17R'};{'17D'};{'ID'}; {'IR'};{'ND'};{'NR'};{'TD'};{'TR'}];%set the names of samples
mkdir('seqs');
%for k = 1:1
        %name = list{k,1};% get the name from the list
        File1 = ['ccs\',name,'.fastq'];% find the file corresponding to the sample
        [header,seqs] = fastqread(File1);% read the headers and sequences from the fastq file
        %% generate a new file containing reads and reversed reads
        [~,len] = size(seqs);
        seqs_c = cell(len*2,2);% create a new cell from headers and sequences
        for i = 1:len%input headers, seqs, and the reverse complimentary seqs to seqs_c
            seq = seqs{1,i};% get one read from the seqs
            seqs_c{i,2} = header{1,i};% set the second column of seqs_c as headers
            seqs_c{len+i,2} = header{1,i};
            seqs_c{i,1}= seq;
            seqs_c{len+i,1}= seqrcomplement(seq);
        end
    save_data(['seqs\',name,'.mat'],seqs_c);
end