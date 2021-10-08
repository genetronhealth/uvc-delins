
all : uvc-rawvcf2hapvcf.out

uvc-rawvcf2hapvcf.out : uvc_rawvcf2hapvcf.cpp
	g++ -std=c++14 -static-libstdc++ -Wall -o uvc-rawvcf2hapvcf.out uvc_rawvcf2hapvcf.cpp /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/ext/htslib-1.11-lowdep/libhts.a -I /bionfsdate/ctDNA/experiment/zhaoxiaofei/uvc/ext/htslib-1.11-lowdep/ -pthread -lm -lz -lbz2 -llzma  -Wextra -g -p

.PHONY clean: 
	rm uvc-rawvcf2hapvcf.out

