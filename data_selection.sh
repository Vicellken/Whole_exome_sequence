# use awk to select data
# select data start with 'rs', then count the number in 3rd column
awk '$3 ~ /^rs/' RawSNP.vcf | awk '{print $3}' | wc -l

# reverse select -- not start with 'rs'
# then select column not start with '#', to remove the annotation data
# count the number in 3rd column
awk '$3 !~ /^rs/' RawSNP.vcf | awk '{if ($1 !~ /^#/) print $3;}' | wc -l