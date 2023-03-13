CXX = g++

CPPFLAGS = -g -std=c++11 -lpthread -ldeflate -lz -O3 -ffast-math -flto=full -DTIMING  -lhts

BUILD_PATH = $(shell pwd)/build

OBJECT = $(BUILD_PATH)/tools.o $(BUILD_PATH)/read.o $(BUILD_PATH)/write.o $(BUILD_PATH)/status.o

LIB_OBJECT = $(BUILD_PATH)/librabbitbamtools.o $(BUILD_PATH)/librabbitbamread.o $(BUILD_PATH)/librabbitbamwrite.o

TARGET = $(BUILD_PATH)/rabbitbam

INCLUDE = /home/user_home/zz/tools/htslib-1.15.1/

LIB = /home/user_home/zz/tools/htslib-1.15.1/build/lib

SHARE = -fPIC -shared




$(BUILD_PATH)/rabbitbam: $(OBJECT) block_mul.cpp
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) block_mul.cpp -o $(BUILD_PATH)/rabbitbam $(OBJECT)

$(BUILD_PATH)/status.o: $(BUILD_PATH)/tools.o Overrepresent.cpp Overrepresent.h  Duplicate.cpp Duplicate.h  BamStatus.cpp BamStatus.h
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB)  Overrepresent.cpp Overrepresent.h  Duplicate.cpp Duplicate.h  BamStatus.cpp BamStatus.h $(SHARE) -o $(BUILD_PATH)/status.o $(BUILD_PATH)/tools.o

$(BUILD_PATH)/write.o : $(BUILD_PATH)/tools.o BamWriter.cpp BamWriter.h BamWriteCompress.h BamWriteCompress.cpp
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamWriter.cpp BamWriter.h BamWriteCompress.h BamWriteCompress.cpp  $(BUILD_PATH)/tools.o $(SHARE) -o $(BUILD_PATH)/write.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamWriter.cpp BamWriter.h BamWriteCompress.h BamWriteCompress.cpp  $(BUILD_PATH)/librabbitbamtools.so $(SHARE) -o $(BUILD_PATH)/librabbitbamwrite.so


$(BUILD_PATH)/read.o : $(BUILD_PATH)/tools.o BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h   BamRead.cpp BamRead.h BamReader.cpp BamReader.h
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h   BamRead.cpp BamRead.h BamReader.cpp BamReader.h $(SHARE) -o $(BUILD_PATH)/read.o  $(BUILD_PATH)/tools.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamCompleteBlock.cpp BamCompleteBlock.h BamCompress.cpp BamCompress.h BamBlock.cpp BamBlock.h   BamRead.cpp BamRead.h BamReader.cpp BamReader.h $(SHARE) -o $(BUILD_PATH)/librabbitbamread.so  $(BUILD_PATH)/librabbitbamtools.so

$(BUILD_PATH)/tools.o :
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamTools.cpp BamTools.h Buffer.cpp Buffer.h config.h $(SHARE) -o $(BUILD_PATH)/tools.o
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB) BamTools.cpp BamTools.h Buffer.cpp Buffer.h config.h $(SHARE) -o $(BUILD_PATH)/librabbitbamtools.so




.PHONY:clean

clean:
	rm -rf $(TARGET) $(OBJECT)
	rm -rf ./output.fq
	rm -rf ./log.txt
	rm -rf ./BamStatus.html
