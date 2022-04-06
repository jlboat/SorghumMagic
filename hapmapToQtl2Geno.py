import sys
import pandas as pd

# will eventually transpose hapmap
with open(sys.argv[1]) as f:
    hapmap = f.read().splitlines()

# ref -> A
# alt -> B
# N -> -

# id,Chr1_81348,Chr1_81600,Chr1_82139,Chr1_83401,Chr1_84018
# indiv1,A,B,A,-,B
# indiv2,A,B,A,-,B

# 11: are samples
# col1=SNP_ID
variants_header = ["id"]
all_variants = []

for x, line in enumerate(hapmap):
    split_line = line.split()
    if x == 0:
        samples = split_line[11:]
    else:
        variants_header.append(split_line[0])
        ref_allele = split_line[1].split('/')[0]
        alt_allele = split_line[1].split('/')[1]
        variant_list = []
        for y, variant in enumerate(split_line[11:]):
            # if y == 0:
            #     variant_list.append(samples[y])
            if variant == ref_allele:
                variant_list.append("A")
            elif variant == alt_allele:
                variant_list.append("B")
            else:
                variant_list.append("-")
        all_variants.append(variant_list)

print(",".join(variants_header))
for x, row in enumerate(pd.DataFrame(all_variants).transpose().values):
    print(samples[x] + "," + ",".join(row))

