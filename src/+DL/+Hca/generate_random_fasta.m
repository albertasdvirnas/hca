function [] = generate_random_fasta(filepath, length)
fid = fopen(filepath, 'w');
fprintf(fid, '>\n');
fprintf('Generating random fasta of length %d bp, please stand by...\n', length)
for i=1:floor(length/70)
  fprintf(fid, randseq(70));
  fprintf(fid, '\n');
end
fprintf(fid, randseq(mod(length, 70)));
fclose(fid);
fprintf('Finished generating random fasta.\n')
end

