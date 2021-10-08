#include <algorithm>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <assert.h>
#include <htslib/faidx.h>
#include <htslib/vcf.h>
#include <htslib/synced_bcf_reader.h>

const auto MIN(const auto a, const auto b) { return ((a) < (b) ? (a) : (b)); }
const auto MAX(const auto a, const auto b) { return ((a) > (b) ? (a) : (b)); }
const auto UPDATE_MIN(auto & a, const auto b) { a = MIN(a, b); }
const auto UPDATE_MAX(auto & a, const auto b) { a = MAX(a, b); }

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

std::vector<std::vector<std::tuple<int, std::string, std::string>>> 
vecof_pos_ref_alt_tup_split(const std::vector<std::tuple<int, std::string, std::string>> & vecof_pos_ref_alt_tup) {
    std::vector<std::vector<std::tuple<int, std::string, std::string>>> vecof_vecof_pos_ref_alt_tup;
    int prev_pos = INT32_MIN;
    for (const auto & pos_ref_alt_tuple : vecof_pos_ref_alt_tup) {
        if (std::get<0>(pos_ref_alt_tuple) >= prev_pos + 6) {
            vecof_vecof_pos_ref_alt_tup.push_back(std::vector<std::tuple<int, std::string, std::string>>());
        }
        prev_pos = std::get<0>(pos_ref_alt_tuple) + (int)MAX(std::get<1>(pos_ref_alt_tuple).size(), std::get<2>(pos_ref_alt_tuple).size());
        vecof_vecof_pos_ref_alt_tup.back().push_back(pos_ref_alt_tuple);
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

int main(int argc, char **argv) {
    const char *fastaref = argv[1];
    const char *uvcvcf = argv[2];
    
    faidx_t *faidx = fai_load(fastaref);
    htsFile *fp = vcf_open(uvcvcf, "r");
    bcf_hdr_t *bcf_hdr = vcf_hdr_read(fp);
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tcHap,Number=1,Type=String,Description=\"Tumor cHap\">");
    bcf_hdr_append(bcf_hdr, "##INFO=<ID=tPRA,Number=1,Type=String,Description=\"Tumor position_REF_ALT, with the three VCF fields separated by underscore\">");
    
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
        char **bcfstring = NULL;
        // float *bcffloats = NULL;
        int32_t *bcfints = NULL;
        
        std::map<std::string, std::vector<std::tuple<int, std::string, std::string>>> map_from_cHap_string_to_vecof_pos_ref_alt_tup;
        
        int nsamples = bcf_hdr_nsamples(bcf_hdr); 
        int sampleidx = nsamples - 1; // last sample
        
        while (bcf_sr_next_line(sr)) {
            bcf1_t *line = bcf_sr_get_line(sr, 0);
            bcf_unpack(line, BCF_UN_ALL);
            ndst_val = 0;
            valsize = bcf_get_format_string(bcf_hdr, line, "cHap", &bcfstring, &ndst_val);
            if (valsize <= 0) { continue; }
            const auto pos_ref_alt_tup = std::make_tuple(line->pos, std::string(line->d.allele[0]), std::string(line->d.allele[1]));
            for (int j = sampleidx; j < nsamples; j++) {
                std::vector<std::string> cHap_substrs = cHapString_to_cHapSubstrs(bcfstring[j]);
                for (const std::string & cHap_string : cHap_substrs) {
                    map_from_cHap_string_to_vecof_pos_ref_alt_tup.insert(std::make_pair(cHap_string, std::vector<std::tuple<int, std::string, std::string>>()));
                    map_from_cHap_string_to_vecof_pos_ref_alt_tup[cHap_string].push_back(pos_ref_alt_tup);
                }
            }
        }
        
        bcf_sr_seek(sr, tname, 0);
        std::cerr << "Will finish processing tname " << tname << "\n";
        std::set<std::tuple<int, std::string, std::string>> complexvar_3tups;
        while (bcf_sr_next_line(sr)) {
            bcf1_t *line = bcf_sr_get_line(sr, 0);
            bcf_unpack(line, BCF_UN_ALL);
            ndst_val = 0;
            valsize = bcf_get_format_string(bcf_hdr, line, "cHap", &bcfstring, &ndst_val);
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
                    if (complexvarDP < 3 || complexvarDP * 5 < vcflineAD) {
                        continue; // this variant has low DP
                    }
                    auto vecof_pos_ref_alt_tup = map_from_cHap_string_to_vecof_pos_ref_alt_tup.find(cHap_string)->second;
                    if (vecof_pos_ref_alt_tup.size() == 0) {
                        continue; // no variant is found
                    }
                    
                    std::sort(vecof_pos_ref_alt_tup.begin(), vecof_pos_ref_alt_tup.end());
                    
                    const auto vecof_vecof_pos_ref_alt_tup = vecof_pos_ref_alt_tup_split(vecof_pos_ref_alt_tup);
                    for (const auto & vecof_pos_ref_alt_tup1 : vecof_vecof_pos_ref_alt_tup) {
                        if (vecof_pos_ref_alt_tup1.size() <= 1) { 
                            continue; // this variant is not complex
                        }
                        
                        int complexvar_begpos = INT32_MAX;
                        int complexvar_endpos = 0;
                        for (auto pos_ref_alt_tup : vecof_pos_ref_alt_tup1) {
                            int pos = std::get<0>(pos_ref_alt_tup);
                            std::string ref = std::get<1>(pos_ref_alt_tup);
                            std::string alt = std::get<2>(pos_ref_alt_tup);
                            int endpos = pos + (int)MAX(alt.size(), ref.size());
                            UPDATE_MIN(complexvar_begpos, pos);
                            UPDATE_MAX(complexvar_endpos, endpos);
                        }
                        std::string complexref = refstring.substr(complexvar_begpos, complexvar_endpos - complexvar_begpos);
                        std::vector<std::string> complexalt_;
                        for (const auto base : complexref) {
                            complexalt_.push_back(std::string(&base, 1));
                        }
                        for (auto pos_ref_alt_tup : vecof_pos_ref_alt_tup1) {
                            int pos = std::get<0>(pos_ref_alt_tup);
                            std::string ref = std::get<1>(pos_ref_alt_tup);
                            std::string alt = std::get<2>(pos_ref_alt_tup);
                            for (int k = pos; k < pos + (int)ref.size(); k++) {
                                complexalt_[k - complexvar_begpos] = "";
                            }
                            complexalt_[pos - complexvar_begpos] = alt;
                        }
                        std::string complexalt = "";
                        for (const auto subalt : complexalt_) {
                            complexalt.append(subalt);
                        }
                        int vcfline_pos = (std::get<0>(pos_ref_alt_tup_from_vcfline));
                        if (complexvar_begpos > vcfline_pos || vcfline_pos >= (complexvar_begpos + (int)MAX(complexref.size(), complexalt.size()))) { 
                            continue; // this var is outside the complex var region
                        }
                        const auto complexvar_3tup = std::make_tuple(complexvar_begpos, complexref, complexalt);
                        
                        if (complexvar_3tups.find(complexvar_3tup) != complexvar_3tups.end()) { 
                            continue; // this var has already been printed
                        }
                        std::cout << tname << "\t" << (complexvar_begpos + 1) << "\t.\t" << complexref << "\t" << complexalt 
                            << "\t.\t.\t" << "tcHap=" << cHap_string 
                            << ";tPRA="
                            << (std::get<0>(pos_ref_alt_tup_from_vcfline)) << "_"
                            << (std::get<1>(pos_ref_alt_tup_from_vcfline)) << "_"
                            << (std::get<2>(pos_ref_alt_tup_from_vcfline))
                            // << "\tGT\t./1"
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
    bcf_hdr_destroy(bcf_hdr);
    vcf_close(fp);
    fai_destroy(faidx);
}

