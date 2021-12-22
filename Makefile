
COMMIT_VERSION=$(shell git rev-parse HEAD | head -c 7)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
COMMIT_DIFF_FULL=$(shell echo "R\"ZXF_specQUOTE(\n $$(git diff HEAD | sed 's/ZXF_specQUOTE/ZXF_specquote/g') \n)ZXF_specQUOTE\"" > gitdiff.txt)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" -DCOMMIT_DIFF_FULL="\"$(COMMIT_DIFF_FULL)\""
CXXFLAGS=-static-libstdc++ -Wall ext/htslib-1.13-lowdep/libhts.a -I ext/htslib-1.13-lowdep/ -pthread -lm -lz -lbz2 -llzma

all : uvcvcf-raw2delins uvcvcf-raw2delins.debug

uvcvcf-raw2delins.debug : uvcvcf_raw2delins.cpp Makefile
	g++ $(VERFLAGS) -O0 -o uvcvcf-raw2delins.debug uvcvcf_raw2delins.cpp $(CXXFLAGS) -p -g -fsanitize=address

uvcvcf-raw2delins : uvcvcf_raw2delins.cpp Makefile
	g++ $(VERFLAGS) -O3 -o uvcvcf-raw2delins uvcvcf_raw2delins.cpp $(CXXFLAGS)

.PHONY: clean deploy

clean:
	rm uvcvcf-raw2delins uvcvcf-raw2delins.debug

deploy:
	cp uvcvcf-raw2delins uvcvcf-raw2delins.debug bin/

