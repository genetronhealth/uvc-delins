
COMMIT_VERSION=$(shell git rev-parse HEAD | head -c 7)
COMMIT_DIFF_SH=$(shell git diff HEAD --shortstat)
COMMIT_DIFF_FULL=$(shell echo "R\"ZXF_specQUOTE(\n $$(git diff HEAD | sed 's/ZXF_specQUOTE/ZXF_specquote/g') \n)ZXF_specQUOTE\"" > gitdiff.txt)
VERFLAGS=-DCOMMIT_VERSION="\"$(COMMIT_VERSION)\"" -DCOMMIT_DIFF_SH="\"$(COMMIT_DIFF_SH)\"" -DCOMMIT_DIFF_FULL="\"$(COMMIT_DIFF_FULL)\""
CXXFLAGS=-static-libstdc++ -Wall ext/htslib-1.13-lowdep/libhts.a -I ext/htslib-1.13-lowdep/ -pthread -lm -lz -lbz2 -llzma

all : uvc-rawvcf2hapvcf.out uvc-rawvcf2hapvcf.debug

uvc-rawvcf2hapvcf.debug : uvc_rawvcf2hapvcf.cpp Makefile
	g++ $(VERFLAGS) -O0 -o uvc-rawvcf2hapvcf.debug uvc_rawvcf2hapvcf.cpp $(CXXFLAGS) -p -g -fsanitize=address

uvc-rawvcf2hapvcf.out : uvc_rawvcf2hapvcf.cpp Makefile
	g++ $(VERFLAGS) -O3 -o uvc-rawvcf2hapvcf.out uvc_rawvcf2hapvcf.cpp $(CXXFLAGS)

.PHONY clean: 
	rm uvc-rawvcf2hapvcf.out uvc-rawvcf2hapvcf.debug

