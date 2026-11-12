awk 'BEGIN {
    FS="[ \t]+"
    OFS="\t"
    print "#FID", "IID", "SEX"
}
NR > 1 && $4 == "male" {
    print $1, $1, 1
}' data/raw/chromosomes/vcf/integrated_call_samples_v3.20130502.ALL.panel \
> data/raw/chromosomes/vcf/1000G_males.psam