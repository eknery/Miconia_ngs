library("seqinr")

f = read.fasta("sequence.fasta")

raw_names = names(f)

names1 = gsub("^.*\\[gene=", "gene=", raw_names)
names2 = gsub("\\].*$", "", names1)
clean_names = gsub("gene=", "", names2)

names(f) = clean_names

write.fasta(f, names = clean_names, file.out ="new.fasta")

clean_names == "rps12"

