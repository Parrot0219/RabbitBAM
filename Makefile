main: 
	g++ -I/home/zz/tools/htslib-1.12/ -L/home/zz/tools/htslib-1.12/build/lib block_mul.cpp BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamRead.cpp BamRead.h BamTools.cpp BamTools.h config.h Overrepresent.cpp Overrepresent.h  Duplicate.cpp Duplicate.h  BamStatus.cpp BamStatus.h  Buffer.cpp Buffer.h BamBlock.cpp BamBlock.h -lpthread -ldeflate -lz  -DTIMING -Ofast -flto  -ffast-math -lhts -o rabbitbam
clear:
	rm -rf ./rabbitbam *.o
clean:
	rm -rf ./rabbitbam *.o
	rm -rf ./output.fq
	rm -rf ./log.txt
	rm -rf ./BamStatus.html
