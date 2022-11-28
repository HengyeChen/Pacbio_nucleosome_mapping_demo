function real = motif_align(aseq,mot_pos,r)

mot = mot_pos.data(r).mot;% read the selected motif
pos = mot_pos.data(r).pos;% read the positions of this motif

size_pos = size(pos);% Check how many motifs on the oligo
len_pos = size_pos(2);
match_rate = zeros(len_pos,1);% create an array for the match rate of motifs
for m = 1:len_pos
    seq = aseq(1, pos(1,m):pos(2,m));% get the sequence at the location of the motif
    [~, align] = swalign(mot,seq,...% align the input sequence to the sequence of motif
        'Alphabet', 'NT',...
        'ScoringMatrix','NUC44');
    match_rate(m) = sum(align(2,:)== '|' | align(2,:)==':')/length(mot);% calculate the match rate at one motif locus
end

sum_rate = sum(match_rate)/len_pos;% calculate the mean match rate of motifs

if sum_rate > 0.99
    real = 1;% if the motifs are perfectly matched, return 1
else
    real = 0;% else return 0
end

