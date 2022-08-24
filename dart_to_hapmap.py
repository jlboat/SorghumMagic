import argparse


def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was designed to convert an Intertek genotype file " +
                                                 "to a HapMap file.\n\n")

    required_named = parser.add_argument_group('required arguments')

    required_named.add_argument(
        "--variants",
        type=str,
        required=True,
        help="The name of the input file (CSV). " +
             "File structure: Intertek genotype format",
        action="store")

    required_named.add_argument(
        "--field-to-geno-ids",
        type=str,
        required=True,
        help="The name of the input file (CSV). " +
             "File structure: Field_ID,Geno_ID",
        action="store")

    return parser.parse_args()


def read_field_id_to_geno(filename):
    field_id_to_geno = {}
    with open(filename) as f:
        for line in f.read().splitlines():
            split_line = line.split(",")
            field_id_to_geno[split_line[0]] = split_line[1]
    return field_id_to_geno


def read_variant_data(filename):
    with open(filename) as f:
        variant_data = f.read().splitlines()
    return variant_data


def get_alleles(line):
    allele_pair = line.split(":")[-1].replace(">", "/")
    return {"alleles": allele_pair,
            "ref": allele_pair.split("/")[0],
            "alt": allele_pair.split("/")[1]}


def unreasonable_function(variant_data, fids_to_gids):
    """A function that should've been split"""
    # S1_1714,T/C,1,1714,+,NA,NA,NA,NA,NA,NA
    first_headers = ["rs", "alleles", "chrom",
                     "pos", "strand", "assembly",
                     "center", "protLSID", "assayLSID",
                     "panelLSID", "QCcode"]
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
                    samples[x] = fids_to_gids[sample]
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
            allele_dict = get_alleles(split_line[0])
            chrom = split_line[chrom_index].replace("Chr0", "").replace("Chr", "")
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
                    string_variants.append(allele_dict["ref"])
                elif value == "1":
                    string_variants.append(allele_dict["alt"])
                else:
                    string_variants.append("N")
            variants = "\t".join(string_variants)
            print("\t".join([snp_name, allele_dict["alleles"], chrom, snp_position, strand, nas, variants]))
    # sys.stderr.write(str(replacements) + "\n")


if __name__ == '__main__':
    arguments = parse_arguments()
    field_ids_to_geno_ids = read_field_id_to_geno(arguments.field_to_geno_ids)
    variant_list = read_variant_data(arguments.variants)
    unreasonable_function(variant_list, field_ids_to_geno_ids)
