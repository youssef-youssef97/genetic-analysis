awk 'BEGIN {
    FS="[ \t]+"
    OFS="\t"
    print "#FID", "IID", "SEX"
}
NR > 1 {
    if ($4 == "male") sex=1
    else if ($4 == "female") sex=2
    else sex="NA"

    print $1, $1, sex
}' data/raw/chromosomes/vcf/integrated_call_samples_v3.20130502.ALL.panel \
> data/raw/chromosomes/vcf/1000G_samples.psam