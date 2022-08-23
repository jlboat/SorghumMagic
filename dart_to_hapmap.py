import sys

seq_id_to_pheno = {}
with open("MAGIC_ID.csv") as f:
    for line in f.read().splitlines():
        split_line = line.split(",")
        seq_id_to_pheno[split_line[0]] = split_line[1]

with open(sys.argv[1]) as f:
    variant_data = f.read().splitlines()

# S1_1714,T/C,1,1714,+,NA,NA,NA,NA,NA,NA 
first_headers = ["rs","alleles","chrom",
        "pos","strand","assembly",
        "center","protLSID","assayLSID",
        "panelLSID","QCcode"]
nas = "NA\tNA\tNA\tNA\tNA\tNA"

replacements = 0
for line in variant_data:
    split_line = line.split(",")
    if line.startswith("*"):
        continue
    elif line.startswith("AlleleID"):
        sample_index = split_line.index("RepAvg") + 1
        samples = split_line[sample_index:]
        for x, sample in enumerate(samples):
            try:
                samples[x] = seq_id_to_pheno[sample]
                replacements += 1
            except KeyError:
                # sys.stderr.write(sample + "\n")
                samples[x] = sample
        header = first_headers + samples
        print("\t".join(header))
        chrom_index = split_line.index("Chrom_Sorghum_v301")
        snp_position_index = split_line.index("ChromPosSnp_Sorghum_v301")
        strand_index = split_line.index("Strand_Sorghum_v301")
    else:
        alleles = split_line[0].split(":")[-1].replace(">","/")
        ref_allele = alleles.split("/")[0]
        alt_allele = alleles.split("/")[1]
        chrom = split_line[chrom_index].replace("Chr0","").replace("Chr","")
        if chrom == "":
            continue
        snp_position = split_line[snp_position_index]
        if snp_position == "0":
            snp_position = split_line[snp_position_index - 1]
        snp_name = "S" + chrom + "_" + snp_position
        strand_string = split_line[strand_index]
        if strand_string.lower() == "minus":
            strand = "-"
        elif strand_string.lower() == "plus":
            strand = "+"
        string_variants = []
        for value in split_line[sample_index:]:
            if value == "0":
                string_variants.append(ref_allele)
            elif value == "1":
                string_variants.append(alt_allele)
            else:
                string_variants.append("N")
        variants = "\t".join(string_variants)
        print("\t".join([snp_name, alleles, chrom, snp_position, strand, nas, variants]))
# sys.stderr.write(str(replacements) + "\n")
