main: 
	g++ -I/home/user_home/zz/tools/htslib-1.15/ -L/home/user_home/zz/tools/htslib-1.15/build/lib block_mul.cpp BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamRead.cpp BamRead.h BamTools.cpp BamTools.h config.h Overrepresent.cpp Overrepresent.h  Duplicate.cpp Duplicate.h  BamStatus.cpp BamStatus.h  Buffer.cpp Buffer.h BamBlock.cpp BamBlock.h -g -std=c++11 -lpthread -ldeflate -lz -O1 -ffast-math -DTIMING  -lhts -o rabbitbam
clear:
	rm -rf ./rabbitbam *.o
clean:
	rm -rf ./rabbitbam *.o
	rm -rf ./output.fq
	rm -rf ./log.txt
	rm -rf ./BamStatus.html
