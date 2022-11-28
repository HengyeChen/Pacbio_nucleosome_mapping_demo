function [C_T_sum,value_sum] = countreads(position,seq_array,value_sum)

[len, wid] = size(seq_array);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        C_T_sum = zeros(len, wid); %C_T_map is an array that contains the sequencing results at each C of GC position    
        for n = 1:len%Read sequence in the cell array one by one
            %convert the format of the sequence from cell to array 
            %celldisp(seq_all);
            C_T_counts = [];
            C_counts = zeros(1, wid);%Count the number of each nucleotide at GC positions
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            for i = 1:length(position(1,:)) %Read the nucleotide at each GC position one by one           
                pos = position(1,i);%Read the posistion(pos) from the position file
                C_T = seq_array(n,pos+1);%Search the nucleotide in one of the sequence
                C_T_counts = [C_T_counts, C_T]; %Output the nucleotide to an array
                %d = 0;
                %1 means C, -1 means T, and 0 means A or T
                if C_T == 'C'
                    C_counts(pos) = 1;
                elseif C_T == 'T'
                    C_counts(pos) = -1;
                    value_sum(2,i) = value_sum(2,i) + C_counts(pos);
                    %d = d+1;
                else
                    C_counts(pos) = 0;          
                end
                
                
            end
            
            C_T_sum(n,:) = C_counts;

        end
        value_sum(2,:) = value_sum(2,:)/len;






