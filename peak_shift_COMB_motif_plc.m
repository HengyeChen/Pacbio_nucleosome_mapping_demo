function [shift_all,shift_dyad,m1] = peak_shift_COMB_motif_plc(seq,mpos,mot_up,mot_down,dyad_f,score_all_pos,ind,bd,c_limit,linker)
%linker = 14;
seq_t = seq;
seq_t(seq_t>0)=0;
seq_t = -seq_t;

shift_old = fix(dyad_f);
shift_dyad = zeros(1,7);
shift_all = zeros(1,1278);
cycle = 0;
while ~isequal(fix(shift_old),fix(shift_dyad)) && cycle <=c_limit
    cycle = cycle + 1;
    shift_all = zeros(1,1278);
    bound = zeros(7,2);
    if sum(shift_dyad)~=0
        shift_old = fix(shift_dyad);
    end
    shift_dyad = zeros(1,7);
    for x = 1:7
        if shift_old(1,x)~=0
            dyad_i = shift_old(1,x);%set the initial dyad position
            %nuc_pass = find(score_all_pos > 0.5);
            bound_down = shift_old(1,x)+bd;
            bound_up = shift_old(1,x)-bd;
            window = 80;
            for u = 1:window %move nucleosome to its upstream
                dyad_up = dyad_i-u;
                if abs(dyad_up-mot_down) > bd || abs(fix(dyad_f(1,x))-mot_down) < bd  %|| ind==144 || ind==145%if a nucleosome is not on motifs at the beginning, do not move it to the top of motifs, unless the sequence is a background sequence(144, 145)
                    if dyad_up-2*bd<=0%if the nucleosome is out of upstream boundary
                        break
                    end
                    if x~=1 && shift_old(1,x-1)~=0 && (dyad_up - shift_old(1,x-1)) < 2*bd + linker
                        break
                    end
                    up_score = score_all_pos(1,dyad_up);
                    if up_score>=score_all_pos(1,dyad_i)%if the new position has same or better score 
                        
                            bound_up = dyad_up-bd;%set upstream boundary
                    else
                        break
                    end
                end
            end
            for d = 1:window%move nucleosome to downstream
                dyad_down = dyad_i+d;
                if abs(dyad_down-mot_up) > bd || abs(fix(dyad_f(1,x))-mot_up) < bd %|| ind==144 || ind==145%if a nucleosome is not on motifs at the beginning, do not move it to the top of motifs, unless the sequence is a background sequence(144, 145)
                    if dyad_down+2*bd>=1278%if the nucleosome is out of downstream boundary
                        break
                    end
                    if x~=7 && shift_old(1,x+1)~=0 && (shift_old(1,x+1)-dyad_down) < 2*bd + linker
                        break
                    end
                    down_score = score_all_pos(1,dyad_down);
                    if down_score>=score_all_pos(1,dyad_i)%if the new position has same or better score
                        bound_down = dyad_down+bd;%set downstream boundary
                    else
                        break
                    end
                end
            end
            %%%%%determine the boundary%%%%%%%%
            if x==1
                bound(x,1)= max(bound_up,1);% boundary position larger than 1
                if shift_old(x+1)~=0
                    bound(x,2)= min(bound_down,(shift_old(x+1)-89));% After move this nucleosome, the distance between it and other nucleosomes needs to be longer than 14bp
                else
                    bound(x,2)= min(bound_down,1200);
                end
            elseif x==7
                bound(x,1)= max(bound_up,shift_old(x-1)+89);
                bound(x,2)= min(bound_down,1278);
            else
                bound(x,1)= max(bound_up,shift_old(x-1)+89);
                if shift_old(1,x+1)~=0
                    bound(x,2)= min(bound_down,(shift_old(x+1)-89));
                else
                    bound(x,2)= min([bound_down,1278]);
                end
            end
            if bound(x,2)-bound(x,1)>2*bd
                bm = mean(bound(x,:));
%                 if bound(x,2)-bound(x,1)>160 && abs(bm - (mpos(1,1)+ mpos(end,end))/2)<=2*bd% if nucleosome larger than 160bp and close to motifs
                    if bm - (mpos(1,1)+ mpos(end,end))/2>0
                        bound(x,1) = bound(x,2)-160;
                    elseif bm - (mpos(1,1)+ mpos(end,end))/2<0
                        bound(x,2) = bound(x,1)+160;
                    end
%                 elseif bound(x,2)-bound(x,1)>160 && abs(bm - (mpos(1,1)+ mpos(end,end))/2)>2*bd% if nucleosome larger than 160bp and far from motifs
%                     bound(x,1) = round(bm)-80;
%                     bound(x,2) = round(bm)+80;
%                 end
                shift_dyad(1,x) = (bound(x,1)+bound(x,2))/2;
                %border = 75;% set the distance between the boundary and the dyad
                %boundary_score = 0.25;
                a = -20; b = (bound(x,2)-bound(x,1))/2;
                if fix(b)==b
                    x1 = a:1:b;
                    m = (a + b)/2;
                    s = 10;
                    pd = makedist('Normal','mu',m,'sigma',s);
                    f = cdf(pd,x1);
                    score = [f(1,22:end), 1, flip(f(1,22:end))];
                else
                    x1 = a:1:round(b);
                    m = (a + b)/2;
                    s = 10;
                    pd = makedist('Normal','mu',m,'sigma',s);
                    f = cdf(pd,x1);
                    score = [f(1,22:end), flip(f(1,22:end))];
                end
                shift_all(1,bound(x,1):bound(x,2)) = (shift_all(1,bound(x,1):bound(x,2))+1).*score;
            else%when nucleosome shorter than 147
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                a = -20; b = bd;
                if fix(b)==b
                    x1 = a:1:b;
                    m = (a + b)/2;
                    s = 10;
                    pd = makedist('Normal','mu',m,'sigma',s);
                    f = cdf(pd,x1);
                    score = [f(1,22:end), 1, flip(f(1,22:end))];
                else
                    x1 = a:1:round(b);
                    m = (a + b)/2;
                    s = 10;
                    pd = makedist('Normal','mu',m,'sigma',s);
                    f = cdf(pd,x1);
                    score = [f(1,22:end), flip(f(1,22:end))];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                shift_dyad(1,x) = shift_old(1,x);
                bound(x,1) = max(shift_dyad(1,x)-bd,1);
                bound(x,2) = min(shift_dyad(1,x)+bd,1278);
                shift_all(1,bound(x,1):bound(x,2)) = (shift_all(1,bound(x,1):bound(x,2))+1).*score;
            end
            
        end
    end
end
m1 = shift_all;
m1(m1<0)=0;












