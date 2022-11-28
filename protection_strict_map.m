function maps = protection_strict_map(seqs)

[len,~]=size(seqs);
maps = zeros(len,1278);
for i = 1:len
    seq = seqs(i,:);
    map = zeros(1,1278);
    GC=find(seq==1);
    C_num= length(GC);
    wid = zeros(1,C_num-1);
    for c = 1:(C_num-1)
        if c == C_num-1
            if GC(1,c) == 1174
            else
                wid(1,c) = 1174-GC(1,c);
                map(1,(GC(1,c)+5):1174) = (zeros(1,(1174-GC(1,c)-5+1))+1);
            end
        else
            if c == 1 && GC(1,1) ~= 27
                wid(1,1) = GC(1,1)-27;
                GT_sub = find(seq(1,1:GC(1,c))==-1);
                if isempty(GT_sub)
                    map(1,(27-5):(27+5)) = (zeros(1,11)+1);
                else
                    map(1,27:(GT_sub(1,end)+5)) = (zeros(1,(GT_sub(1,end)+5-27+1))+1);
                end
            end
            wid(1,c) = GC(1,c+1)-GC(1,c);
            %         if wid(1,c)>30
            GT_sub = find(seq(1,GC(1,c):GC(1,c+1))==-1);
            if isempty(GT_sub)~=1
                %sscore = smooth_peak(-20,5);
                map(1,(GC(1,c)-5):(GC(1,c)+5))= zeros(1,11)-1;%).*sscore;
                
                if length(GT_sub) == 1
                    map(1,(GC(1,c)+(GT_sub(1,1)-5)):(GC(1,c)+(GT_sub(1,1)+5)))= map(1,(GC(1,c)+(GT_sub(1,1)-5)):(GC(1,c)+(GT_sub(1,1)+5)))+1;
                else
                    %tscore = smooth_peak(-20,((GC(1,c)+(GT_sub(1,end)+5))-max([(GC(1,c)+(GT_sub(1,1)-5)),1]))/2);
                    map(1,max([(GC(1,c)+(GT_sub(1,1)-5)),1]):(GC(1,c)+(GT_sub(1,end)+5)))= (map(1,max([(GC(1,c)+(GT_sub(1,1)-5)),1]):(GC(1,c)+(GT_sub(1,end)+5)))+1);
                end
            else
                %c_bound = [GC(1,c),GC(1,c+1)];
                %cscore = smooth_peak(-20,(GC(1,c+1)-GC(1,c))/2);
                map(1,GC(1,c):GC(1,c+1)) = (zeros(1,(GC(1,c+1)-GC(1,c)+1))-1);
            end
        end
    end
    maps(i,:) = map;
end
i=1;







