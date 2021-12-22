#!/usr/bin/env sh

if [ "${1}" = '-h' -o $# -lt 3 ]; then
    printf "Usage: ${0} <REF-FASTA> <INPUT-VCF-GZ> <OUTPUT-PREFIX> <OTHER-PARAMS...> \n"
    printf "  REF-FASTA: reference fasta file \n"
    printf "  INPUT-VCF-GZ: input vcf file in block gzip format\n"
    printf "  OUTPUT-PREFIX: the prefix of the output VCF files. Each output VCF file has the following suffix.\n"
    printf "    .consumed-simple-with-fmt.vcf.gz : vcf containing simple variants that should be removed because they are entirely parts of some complex variants. \n"
    printf "    .leftover-simple-with-fmt.vcf.gz : vcf containing simple variants that are not consumed by any complex variant. \n"
    printf "    .generated-delins.vcf.gz         : vcf containing delins variants combined from simple variants. \n"
    printf "    .leftover-simple.vcf.gz          : same as unhapvcf but with sample format data removed. \n"
    printf "    .merged-simple-delins.vcf.gz     : vcf containing both simple and delins variants. \n"
    exit -1
fi

scriptdir="$(dirname $(which "${0}"))"
UVCdelins="${scriptdir}/uvcvcf-raw2delins"
FASTAREF="${1}"
rawvcf="${2}"
combvcf="${3}.consumed-simple-with-fmt.vcf.gz" # vcf containing simple variants that should be removed because they are entirely parts of some complex variants.
unhapvcf="${3}.leftover-simple-with-fmt.vcf.gz" # vcf containing simple variants that are not consumed by any complex variant
hapvcf="${3}.generated-delins.vcf.gz" # vcf containing delins variants combined from simple variants
unhapvcf2="${3}.leftover-simple.vcf.gz" # same as unhapvcf but with sample format data removed
mergedvcf="${3}.merged-simple-delins.vcf.gz" # vcf containing both simple and delins variants

bcftools index -f "${rawvcf}" || true

"${UVCdelins}" "${FASTAREF}" "${rawvcf}" -C "${combvcf}" -D "${unhapvcf}" -M wz "${@:4}" 2> "${rawvcf/uvc.vcf.gz/uvc-hap.stderr}" \
     |   bcftools sort - \
     |   bcftools view -Oz -o "${hapvcf}" -

bcftools index -f "${combvcf}"  || true
bcftools index -f "${unhapvcf}" || true
bcftools index -f "${hapvcf}"   || true
bcftools view -s '' --force-samples -Oz -o "${unhapvcf2}" "${unhapvcf}"
bcftools index -f "${unhapvcf2}" || true
bcftools merge -Oz -m none -o "${mergedvcf}" "${unhapvcf2}" "${hapvcf}"
bcftools index -f "${mergedvcf}" || true # step-uvc-hap-1

