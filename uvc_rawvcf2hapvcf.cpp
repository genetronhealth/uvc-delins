#include <algorithm>
#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <assert.h>
#include <float.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

#include <math.h>
#include <unistd.h>

const auto MIN(const auto a, const auto b) { return ((a) < (b) ? (a) : (b)); }
const auto MAX(const auto a, const auto b) { return ((a) > (b) ? (a) : (b)); }
const auto UPDATE_MIN(auto & a, const auto b) { a = MIN(a, b); }
const auto UPDATE_MAX(auto & a, const auto b) { a = MAX(a, b); }

class VariantInfo {
public:
    float qual;
    int tbDP;
    int tDP;
    std::array<int, 2> tADR;
    VariantInfo(float q, int dp1, int dp2, std::array<int, 2> tADR1) {
        qual = q;
        tbDP = dp1;
        tDP = dp2;
        for (int i = 0; i < 2; i++) {
            tADR[i] = tADR1[i];
        }
    };
    bool operator < (const VariantInfo & vi2) const {
        return 0;
    }
};

template <class T>
std::string
other_join(const T & container, std::string sep = std::string(",")) {
    std::string ret = "";
    for (const auto & e : container) {
        ret += std::to_string(e) + sep;
    }
    if (ret.size() > 0) { ret.pop_back(); }
    return ret;
}

int 
cHapSubstr_to_totDP(const std::string & cHapSubstr) {
    int ret = 0;
    int tot_n_amperstands = 0;
    for (const auto ch : cHapSubstr) {
        if ('&' == ch) {
            tot_n_amperstands++;
        }
    }
    std::string tmp;
    std::stringstream ss(cHapSubstr);
    int ntokens = 0;
    while(getline(ss, tmp, '&')) {
        if ((tot_n_amperstands - 1 == ntokens) || (tot_n_amperstands == ntokens)) {
            ret += atoi(tmp.c_str());
        }
        ntokens++;
    }
    return ret;
}

inline int varlen2reflen(int varlen) {
    if (varlen > 0) { return varlen * 2; }
    else { return -varlen; }
}

std::vector<std::vector<std::tuple<int, std::string, std::string, VariantInfo>>> 
vecof_pos_ref_alt_tup_split(const std::vector<std::tuple<int, std::string, std::string, int, int, VariantInfo>> & vecof_pos_ref_alt_begpos_endpos_tup,
        int defaultCB, int defaultCO, int defaultCE, bool enable_short_tandem_repeat_adjust) {
    std::vector<std::vector<std::tuple<int, std::string, std::string, VariantInfo>>> vecof_vecof_pos_ref_alt_tup;
    int prev_pos = INT32_MIN;
    int prev_varlen = 0;
    int nextof_prev = 0;
    int prevof_curr = 0;
    int prev_link_n_bases = 0;
    int delim_pos = 0;
    for (const auto & pos_ref_alt_begpos_endpos_tuple : vecof_pos_ref_alt_begpos_endpos_tup) {
        int varlen = ((int)(std::get<1>(pos_ref_alt_begpos_endpos_tuple).size()) - (int)(std::get<2>(pos_ref_alt_begpos_endpos_tuple).size()));
        const VariantInfo & variantinfo = std::get<5>(pos_ref_alt_begpos_endpos_tuple);
        
        int link_n_bases = 0;
        if (0 == varlen && 0 == prev_varlen) {
            link_n_bases = defaultCB;
        } else {
            link_n_bases = defaultCO + (MAX(varlen2reflen(varlen), abs(prev_varlen)) * defaultCE);
        }
        UPDATE_MAX(delim_pos, prev_pos + MAX(link_n_bases, prev_link_n_bases));
        int curr_pos = std::get<0>(pos_ref_alt_begpos_endpos_tuple);
        if (curr_pos >= delim_pos && ((!enable_short_tandem_repeat_adjust) || nextof_prev <= prevof_curr)) {
            vecof_vecof_pos_ref_alt_tup.push_back(std::vector<std::tuple<int, std::string, std::string, VariantInfo>>());
        }
        prev_link_n_bases = link_n_bases;
        prev_pos = std::get<0>(pos_ref_alt_begpos_endpos_tuple) + (int)MAX(std::get<1>(pos_ref_alt_begpos_endpos_tuple).size(), std::get<2>(pos_ref_alt_begpos_endpos_tuple).size());
        prev_varlen = varlen;
        nextof_prev = std::get<4>(pos_ref_alt_begpos_endpos_tuple);
        vecof_vecof_pos_ref_alt_tup.back().push_back(
            std::make_tuple(
                std::get<0>(pos_ref_alt_begpos_endpos_tuple),
                std::get<1>(pos_ref_alt_begpos_endpos_tuple),
                std::get<2>(pos_ref_alt_begpos_endpos_tuple),
                variantinfo
            )
        );
    }
    return vecof_vecof_pos_ref_alt_tup;
}

const std::vector<std::string> 
cHapString_to_cHapSubstrs(const std::string & cHapString) {
    std::vector<std::string> cHap_substrs;
    std::string cHap_substr = "";
    int parenlevel = 0;
    for (const auto cHap_char : cHapString) {
        if ('(' == cHap_char) {
            parenlevel++; 
        }
        if (')' == cHap_char) {
            parenlevel--; 
        }
        cHap_substr.push_back(cHap_char);
        if (0 == parenlevel) {
            cHap_substrs.push_back(cHap_substr);
            cHap_substr = "";
        }
    }
    return cHap_substrs;
}

std::map<std::string, int> build_tname2tid_from_faidx(const faidx_t *faidx) {
    std::map<std::string, int> ret;
    for (int i = 0; i < faidx_nseq(faidx); i++) {
        ret.insert(std::make_pair(faidx_iseq(faidx, i), i));
    }
    return ret;
}

const int DEFAULT_D = 3;
const double DEFAULT_F = 0.1 + 1e-6;
const std::string DEFAULT_H = "cHap";
const int DEFAULT_CB = 4; // SNV to SNV
const int DEFAULT_CO = 6; // (SNV to InDel gap-open) and (InDel to InDel gap-open)
const int DEFAULT_CE = 1; // (SNV to InDel gap-ext ) and (InDel to InDel gap-ext )
// const int DEFAULT_CR = 10;   // (SNV to InDel repeat-decrease) and (InDel to InDel repeat-decrease)
const double  POWLAW_EXPONENT = 3.0;

void help(int argc, char **argv) {
    fprintf(stderr, "Program %s version %s ( %s )\n", argv[0], COMMIT_VERSION, COMMIT_DIFF_SH);
    fprintf(stderr, "  This program combines simple variants into complex variants. \n");
    
    fprintf(stderr, "Usage: %s <REFERENCE-FASTA> <UVC-VCF-GZ> \n", argv[0]);
    fprintf(stderr, "Optional parameters:\n");
    fprintf(stderr, " -d minimum allele depth of the linked variants [default to %d].\n", DEFAULT_D);
    fprintf(stderr, " -f minimum fraction of the linked variants [default to %f].\n", DEFAULT_F);
    fprintf(stderr, " -h FORMAT tag in the UVC-VCF-GZ file used to contain the haplotype information [default to %s].\n", DEFAULT_H.c_str());
    fprintf(stderr, " -p the power-law exponent for computing tumor haplotype variant qualities [default to %f].\n", POWLAW_EXPONENT);    
    fprintf(stderr, " -B maximum number of bases between SNV and SNV to be considered as linked [default to %d].\n", DEFAULT_CB);
    fprintf(stderr, " -O gap opening for the maximum number of bases between InDel and SNV/InDel to be considered as linked [default to %d].\n", DEFAULT_CO);
    fprintf(stderr, " -E gap extension for the maximum number of bases between InDel and SNV/InDel to be considered as linked [default to %d].\n", DEFAULT_CE);
    fprintf(stderr, " -T bed file that overrides the -d, -f, -B -O, and -E parameters in the defined regions [default to None].\n");
    
    fprintf(stderr, " -S boolean flag indicating if short-tandem-repeats (STRs) should be considered in the merging of simple variants [default to false].\n");
    fprintf(stderr, " -L boolean flag indicating if left-trimming of bases occurring in both REF and ALT should be disabled [default to false].\n");
    fprintf(stderr, " -R boolean flag indicating if right-trimming of bases occurring in both REF and ALT should be disabled [default to false].\n");
    
    exit(-1);
}

int main(int argc, char **argv) {
    
    char *fastaref = NULL;
    char *uvcvcf = NULL;
    char *bedfile = NULL;
    double powlaw_exponent = POWLAW_EXPONENT;
    int linkdepth1 = DEFAULT_D;
    double linkfrac1 = DEFAULT_F;
    std::string defaultH1 = DEFAULT_H;
    int defaultCB1 = DEFAULT_CB;
    int defaultCO1 = DEFAULT_CO;
    int defaultCE1 = DEFAULT_CE;
    bool enable_short_tandem_repeat_adjust = false;
    bool disable_left_trim = false;
    bool disable_right_trim = false;
    int opt = -1;
    while ((opt = getopt(argc, argv, "b:d:f:p:B:O:E:S:T:LR")) != -1) {
        switch (opt) {
            case 'd': linkdepth1 = atoi(optarg); break;
            case 'f': linkfrac1 = atof(optarg); break;
            case 'h': defaultH1 = optarg; break;
            case 'p': powlaw_exponent = atof(optarg); break;
            case 'T': bedfile = optarg; break;
            case 'B': defaultCB1 = atoi(optarg); break;
            case 'O': defaultCO1 = atoi(optarg); break;
            case 'E': defaultCE1 = atoi(optarg); break;
            case 'S': enable_short_tandem_repeat_adjust = true; break;
            case 'L': disable_left_trim = true; break;
            case 'R': disable_right_trim = true; break;
            default: help(argc, argv);
        }
    }
    for (int posidx = 0; optind < argc; optind++, posidx++) {
        if (0 == posidx) { fastaref = argv[optind]; }
        else if (1 == posidx) { uvcvcf = argv[optind]; }
    }
    if (NULL == fastaref || NULL == uvcvcf) {
        help(argc, argv);
    }
    
    faidx_t *faidx = fai_load(fastaref);
    htsFile *fp = vcf_open(uvcvcf, "r");
    bcf_hdr_t *bcf_hdr = vcf_hdr_read(fp);
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tHap,Number=1,Type=String,Description=\"Tumor [x]Hap where [x] can be b, c, or c2.\">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tPRA,Number=1,Type=String,Description=\"Tumor position_REF_ALT, with the three VCF fields separated by underscore\">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tDPm,Number=1,Type=Integer,Description=\"Tumor total deduped depth (deprecated, please see CDP1f and CDP1r) taken as the minimum of the decomposed SNV-InDel alleles. \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tDPM,Number=1,Type=Integer,Description=\"Tumor total deduped depth (deprecated, please see CDP1f and CDP1r) taken as the maximum of the decomposed SNV-InDel alleles. \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tADA,Number=A,Type=Integer,Description=\"Tumor total deduped depth of each MNV or complex variant. \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tADRm,Number=R,Type=Integer,Description=\"Tumor deduped depth of each MNV or complex variant by using the minimum depth among the linked SNVs and/or InDels. \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tADRM,Number=R,Type=Integer,Description=\"Tumor deduped depth of each MNV or complex variant by using the maximum depth among the linked SNVs and/or InDels (deprecated, please see cDP1f and cDP1r) (WARNING: use this field with caution because it should not be used under normal circumstances!). \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tAD2F,Number=A,Type=Integer,Description=\"Percentage (100x) of reads that support the complex variant among the decomposed non-complex variants. \">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tHVQ,Number=A,Type=Integer,Description=\"Haplotype (non-SNV and non-InDel small) variant Quality. \">");
    
    bcf_hdr_t *bcf_hdr2 = bcf_hdr_dup(bcf_hdr);
    // int set_samples_ret = bcf_hdr_set_samples(bcf_hdr, bcf_hdr->samples[bcf_hdr->nsamples_ori - 1], false);
    int set_samples_ret1 = bcf_hdr_set_samples(bcf_hdr2, NULL, false);
    assert(0 == set_samples_ret1);
    kstring_t hdr_str = {0, 0, NULL};
    bcf_hdr_format(bcf_hdr2, true, &hdr_str);
    std::cout << hdr_str.s;
    bcf_hdr_destroy(bcf_hdr2);
    // int set_samples_ret2 = bcf_hdr_set_samples(bcf_hdr, "-", false);
    // assert(0 == set_samples_ret2);
    
    int vcf_nseqs = -1;
    const char **seqnames = bcf_hdr_seqnames(bcf_hdr, &vcf_nseqs);
    
    assert (faidx_nseq(faidx) == vcf_nseqs);
    
    for (int i = 0; i < faidx_nseq(faidx); i++) {
        const char *seqname = faidx_iseq(faidx, i);
        assert(!strcmp(seqname, seqnames[i])); 
    }
    
    std::ifstream bedstream;
    std::map<std::string, int> tname2tid;
    if (NULL != bedfile) { 
        bedstream.open(bedfile); 
        tname2tid = build_tname2tid_from_faidx(faidx);
    }
    int bedtid = -1;
    int bedbeg = -1;
    int bedend = -1;

    for (int i = 0; i < vcf_nseqs; i++) {
        const char *tname = faidx_iseq(faidx, i);
        int regionlen;
        char *fetchedseq = fai_fetch(faidx, tname, &regionlen);
        const std::string refstring = fetchedseq;
        free(fetchedseq);
        
        bcf_srs_t *const sr = bcf_sr_init();
        if (NULL == sr) {
            std::cerr << "Failed to initialize bcf sr\n";
            exit(-6);
        }

        bcf_sr_set_regions(sr, tname, false);
        int sr_set_opt_retval = bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
        if (sr_set_opt_retval < 0) {
            exit(-8);
        }
        int sr_add_reader_retval = bcf_sr_add_reader(sr, uvcvcf);
        if (sr_add_reader_retval != 1) {
            exit(-9);
        }

        int valsize = 0;
        int ndst_val = 0;
        int string_valsize = 0;
        int string_ndst_val = 0;
        char **bcfstring = NULL;
        int32_t *bcfints = NULL;
        
        std::map<std::string, std::vector<std::tuple<int, std::string, std::string, int, int, VariantInfo>>> map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup;
        
        int nsamples = bcf_hdr_nsamples(bcf_hdr); 
        int sampleidx = nsamples - 1; // last sample
        
        while (bcf_sr_next_line(sr)) {
            bcf1_t *line = bcf_sr_get_line(sr, 0);
            bcf_unpack(line, BCF_UN_ALL);
            ndst_val = 0;
            
            fprintf(stderr, "Processing the line: tid = %d pos = %ld\n", line->rid, line->pos);
            valsize = bcf_get_info_int32(bcf_hdr, line, "R3X2", &bcfints, &ndst_val);
            if (valsize < 6) { continue; }
            const int posleft = bcfints[0];
            const int posright = bcfints[3] + (bcfints[4] * bcfints[5]);
            
            if (bcfstring) {
                free(bcfstring[0]);
                free(bcfstring);
                bcfstring = NULL;
                string_ndst_val = 0;
            }
            string_valsize = bcf_get_format_string(bcf_hdr, line, defaultH1.c_str(), &bcfstring, &string_ndst_val);
            assert (1 <= string_valsize || !fprintf(stderr, "The size of cHap is %d instead of at least 1!\n", string_valsize));
            valsize = bcf_get_info_int32(bcf_hdr, line, "tbDP", &bcfints, &ndst_val);
            int tbDP = bcfints[0];
            valsize = bcf_get_info_int32(bcf_hdr, line, "tDP", &bcfints, &ndst_val);
            int tDP = bcfints[0];
            valsize = bcf_get_info_int32(bcf_hdr, line, "tADR", &bcfints, &ndst_val);
            std::array<int, 2> tADR = std::array<int, 2>({bcfints[0], bcfints[1]});
            
            const auto pos_ref_alt_begpos_endpos_tup = std::make_tuple(line->pos, std::string(line->d.allele[0]), std::string(line->d.allele[1]), posleft, posright,
                    VariantInfo(line->qual, tbDP, tDP, tADR));
            assert (tDP > 0 || !fprintf(stderr, "%d > 0 failed for rid - %d pos - %ld ref - %s alt - %s!\n", tDP, line->rid, line->pos, line->d.allele[0], line->d.allele[1]));
            for (int j = sampleidx; j < nsamples; j++) {
                std::vector<std::string> cHap_substrs = cHapString_to_cHapSubstrs(bcfstring[j]);
                for (const std::string & cHap_string : cHap_substrs) {
                    map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup.insert(std::make_pair(cHap_string, std::vector<std::tuple<int, std::string, std::string, int, int, VariantInfo>>()));
                    map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup[cHap_string].push_back(pos_ref_alt_begpos_endpos_tup);
                }
            }
        }
        
        bcf_sr_seek(sr, tname, 0);
        std::cerr << "Will finish processing tname " << tname << "\n";
        std::set<std::tuple<int, std::string, std::string>> complexvar_3tups;
        while (bcf_sr_next_line(sr)) {
           
            int linkdepth = linkdepth1;
            int linkfrac = linkfrac1;
            int defaultCB = defaultCB1;
            int defaultCO = defaultCO1;
            int defaultCE = defaultCE1;
            
            bcf1_t *line = bcf_sr_get_line(sr, 0);
            bcf_unpack(line, BCF_UN_ALL); 

            if (bedfile != NULL && bedstream.good()) {
                int vcftid = line->rid;
                int vcfbeg = line->pos;
                while (bedtid < vcftid || (bedtid == vcftid && bedbeg < vcfbeg)) {
                    std::string line;
                    getline(bedstream, line);
                    if (line.empty()) { break; }
                    std::istringstream linestream(line);
                    std::string bedchrom;
                    linestream >> bedchrom;
                    linestream >> bedbeg;
                    linestream >> bedend;
                    if (tname2tid.find(bedchrom) == tname2tid.end()) {
                        fprintf(stderr, "The bed file %s with line %s has invalid tname %s so this line is skipped!\n", bedfile, line.c_str(), bedchrom.c_str());
                        bedtid = -1;
                    } else {
                        bedtid = tname2tid[bedchrom];
                        std::string token;
                        while (linestream.good()) {
                            linestream >> token;
                            if (!token.compare("-d")) {
                                linestream >> token;
                                linkdepth = atoi(token.c_str());
                            }
                            if (!token.compare("-f")) {
                                linestream >> token;
                                linkfrac = atof(token.c_str());
                            }
                            if (!token.compare("-B")) {
                                linestream >> token;
                                defaultCB = atoi(token.c_str());
                            }
                            if (!token.compare("-O")) {
                                linestream >> token;
                                defaultCO = atoi(token.c_str());
                            }
                            if (!token.compare("-E")) {
                                linestream >> token;
                                defaultCE = atoi(token.c_str());
                            }
                        }
                    }
                }
            }
            
            ndst_val = 0;
            valsize = bcf_get_format_string(bcf_hdr, line, defaultH1.c_str(), &bcfstring, &ndst_val);
            if (valsize <= 0) { continue; }
            ndst_val = 0;
            valsize = bcf_get_format_int32(bcf_hdr, line, "AD", &bcfints, &ndst_val);
            if (valsize <= 0) { continue; }
            const int vcflineAD = bcfints[ndst_val - 1];
            const auto pos_ref_alt_tup_from_vcfline = std::make_tuple(line->pos, std::string(line->d.allele[0]), std::string(line->d.allele[1]));
            
            for (int j = sampleidx; j < nsamples; j++) {
                std::vector<std::string> cHap_substrs = cHapString_to_cHapSubstrs(bcfstring[j]);
                for (std::string cHap_string : cHap_substrs) {
                    int complexvarDP = cHapSubstr_to_totDP(cHap_string);
                    if (complexvarDP < linkdepth || complexvarDP < linkfrac * vcflineAD) {
                        // fprintf(stderr, "Skipping the variant %s %d because it has low DP\n", tname, line->pos);
                        continue; // this variant has low DP
                    }
                    auto vecof_pos_ref_alt_begpos_endpos_tup = map_from_cHap_string_to_vecof_pos_ref_alt_begpos_endpos_tup.find(cHap_string)->second;
                    if (vecof_pos_ref_alt_begpos_endpos_tup.size() == 0) {
                        // fprintf(stderr, "Skipping the variant %s %d because it has no variants!\n", tname, line->pos);
                        continue; // no variant is found
                    }
                    
                    std::sort(vecof_pos_ref_alt_begpos_endpos_tup.begin(), vecof_pos_ref_alt_begpos_endpos_tup.end());
                    
                    const auto vecof_vecof_pos_ref_alt_tup = vecof_pos_ref_alt_tup_split(vecof_pos_ref_alt_begpos_endpos_tup,
                        defaultCB, defaultCO, defaultCE, enable_short_tandem_repeat_adjust);
                    for (const auto & vecof_pos_ref_alt_tup1 : vecof_vecof_pos_ref_alt_tup) {
                        
                        if (vecof_pos_ref_alt_tup1.size() <= 1) { 
                            // fprintf(stderr, "Skipping the variant %s %d because it is not complex\n", tname, line->pos);
                            continue; // this variant is not complex
                        }
                        double errmodel_complexvarfrac_phred = 0.0;
                        int complexvar_begpos = INT32_MAX;
                        int complexvar_endpos = 0;
                        float cv_qual = FLT_MAX;
                        int tDPmin = INT_MAX;
                        int tDPmax = 0;
                        std::array<int, 2> tADRmin = {INT_MAX, INT_MAX};
                        std::array<int, 2> tADRmax = {0};
                        for (auto pos_ref_alt_tup : vecof_pos_ref_alt_tup1) {
                            int pos = std::get<0>(pos_ref_alt_tup);
                            std::string ref = std::get<1>(pos_ref_alt_tup);
                            std::string alt = std::get<2>(pos_ref_alt_tup);
                            int endpos = pos + (int)MAX(alt.size(), ref.size());
                            UPDATE_MIN(complexvar_begpos, pos);
                            UPDATE_MAX(complexvar_endpos, endpos);
                            float qual = std::get<3>(pos_ref_alt_tup).qual;
                            int tDP = std::get<3>(pos_ref_alt_tup).tDP;
                            const auto &tADR = std::get<3>(pos_ref_alt_tup).tADR;
                            
                            UPDATE_MIN(tDPmin, tDP);
                            UPDATE_MAX(tDPmax, tDP);
                            for (int j = 0; j < 2; j++) {
                                UPDATE_MIN(tADRmin[j], tADR[j]);
                                UPDATE_MAX(tADRmax[j], tADR[j]);
                                assert (tADR[1] <= tDP);
                                assert (complexvarDP <= tADR[1]);
                            }
                            errmodel_complexvarfrac_phred += 10.0 / log(10.0) * log((double)(tADR[1] - complexvarDP + 0.5) / (double)(tDP + 1.0));
                            cv_qual = MIN(qual, cv_qual);
                        }
                        double complexvarfrac = (double)(complexvarDP + 0.5) / (double)(tDPmax + 1.0);
                        double complexvarfrac_ratio_phred = 10.0 / log(10.0) * log(complexvarfrac) - errmodel_complexvarfrac_phred;
                        
                        std::string complex_ref = refstring.substr(complexvar_begpos, complexvar_endpos - complexvar_begpos);
                        std::vector<std::string> complex_alt_;
                        for (const auto base : complex_ref) {
                            complex_alt_.push_back(std::string(&base, 1));
                        }
                        for (auto pos_ref_alt_tup : vecof_pos_ref_alt_tup1) {
                            int pos = std::get<0>(pos_ref_alt_tup);
                            std::string ref = std::get<1>(pos_ref_alt_tup);
                            std::string alt = std::get<2>(pos_ref_alt_tup);
                            for (int k = pos; k < pos + (int)ref.size(); k++) {
                                complex_alt_[k - complexvar_begpos] = "";
                            }
                            complex_alt_[pos - complexvar_begpos] = alt;
                        }
                        std::string complex_alt = "";
                        for (const auto subalt : complex_alt_) {
                            complex_alt.append(subalt);
                        }
                        int vcfline_pos = (std::get<0>(pos_ref_alt_tup_from_vcfline));
                        if (complexvar_begpos > vcfline_pos || vcfline_pos >= (complexvar_begpos + (int)MAX(complex_ref.size(), complex_alt.size()))) { 
                            /*
                            fprintf(stderr, "Skipping the variant %s %d %s %s because it is outside the ROI\n", 
                                tname, vcfline_pos, 
                                std::get<1>(pos_ref_alt_tup_from_vcfline).c_str(), 
                                std::get<2>(pos_ref_alt_tup_from_vcfline).c_str());
                            */
                            continue; // this var is outside the complex var region
                        }
                        const auto complexvar_3tup = std::make_tuple(complexvar_begpos, complex_ref, complex_alt);
                        
                        if (complexvar_3tups.find(complexvar_3tup) != complexvar_3tups.end()) { 
                            fprintf(stderr, "Skipping the variant %s %d %s %s because it is already visited\n", 
                                tname, vcfline_pos, 
                                std::get<1>(pos_ref_alt_tup_from_vcfline).c_str(), 
                                std::get<2>(pos_ref_alt_tup_from_vcfline).c_str());
                            continue; // this var has already been printed
                        }
                        int cv_begpos = complexvar_begpos;
                        std::string cv_ref = complex_ref;
                        std::string cv_alt = complex_alt;
                        if (!disable_left_trim) {
                            size_t begpos_inc = 0;
                            while ((begpos_inc + 1 < cv_ref.size())
                                    && (begpos_inc + 1 < cv_alt.size())
                                    && (cv_ref[begpos_inc] == (cv_alt[begpos_inc]))) {
                                begpos_inc++;
                            }
                            cv_begpos += begpos_inc;
                            cv_ref = cv_ref.substr(begpos_inc);
                            cv_alt = cv_alt.substr(begpos_inc);
                        }
                        if (!disable_right_trim) {
                            size_t endpos_dec = 0;
                            while (cv_ref.size() > (1 + endpos_dec) && cv_alt.size() > (1 + endpos_dec)
                                    && cv_ref[cv_ref.size() - (1 + endpos_dec)] == cv_alt[cv_alt.size() - (1 + endpos_dec)]) {
                                endpos_dec++;
                            }
                            cv_ref = cv_ref.substr(0, cv_ref.size() - endpos_dec);
                            cv_alt = cv_alt.substr(0, cv_alt.size() - endpos_dec);
                        }
                        std::cout << tname << "\t" << (cv_begpos + 1) << "\t.\t" << cv_ref << "\t" << cv_alt 
                            << "\t" << std::to_string(cv_qual) << "\t.\t" << "tHap=" << cHap_string 
                            << ";tPRA="
                            << (std::get<0>(pos_ref_alt_tup_from_vcfline)) << "_"
                            << (std::get<1>(pos_ref_alt_tup_from_vcfline)) << "_"
                            << (std::get<2>(pos_ref_alt_tup_from_vcfline))
                            << ";tDPm=" << tDPmin
                            << ";tDPM=" << tDPmax
                            << ";tADA=" << complexvarDP
                            << ";tADRm=" << other_join(tADRmin, ",")
                            << ";tADRM=" << other_join(tADRmax, ",")
                            << ";tAD2F=" << (tADRmin[1] * 100 / MAX(1, tADRmax[1]))
                            << ";tHVQ=" << (powlaw_exponent * complexvarfrac_ratio_phred)
                            << "\n";
                        for (auto it = complexvar_3tups.begin(); it != complexvar_3tups.end(); ) {
                            int endpos = std::get<0>(*it) + (int)MAX(std::get<1>(*it).size(), std::get<2>(*it).size());
                            if (endpos < vcfline_pos) {
                                it = complexvar_3tups.erase(it);
                            } else {
                                it++;
                            }
                        }
                        complexvar_3tups.insert(complexvar_3tup);
                    }
                }
            }
        }
        bcf_sr_destroy(sr);
    }
    if (NULL != bedfile) { 
        bedstream.close(); 
    }
    bcf_hdr_destroy(bcf_hdr);
    vcf_close(fp);
    fai_destroy(faidx);
}
