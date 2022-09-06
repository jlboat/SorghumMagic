# SorghumMagic
Scripts and methods used to analyze the Sorghum MAGIC population

```
.
├── data
│   ├── dart
│   │   ├── Intertek.hapmap.txt
│   │   ├── Intertek.snps.csv
│   │   └── MAGIC_ID.csv
│   ├── phenotypes.csv
│   └── qtl2
├── LICENSE
├── README.md
└── scripts
    ├── convert
    │   └── dart_to_hapmap.py
    ├── gwas
    │   └── GWAS_MAGIC.R
    ├── misc
    │   └── Correlation_Heatplot.R
    └── qtl_mapping
        ├── estimate_geneticmap.R
        └── qtl2_analysis.R
```

## File conversion

```
## Convert Intertek file to HapMap (also converts field IDs to genotype IDs
python ./scripts/convert/dart_to_hapmap.py \
    --variants ./data/dart/Intertek.snps.csv \
    --field-to-geno-ids ./data/dart/MAGIC_ID.csv > ./data/dart/Intertek.hapmap.txt 
```
