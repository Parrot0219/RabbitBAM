CXX = g++

CPPFLAGS = -g -std=c++11 -lpthread -ldeflate -lz -O3 -ffast-math -flto=full -DTIMING  -lhts

BUILD_PATH = /home/user_home/zz/API/BamUniversalStatus/build

OBJECT = $(BUILD_PATH)/tools.o $(BUILD_PATH)/read.o $(BUILD_PATH)/write.o $(BUILD_PATH)/status.o

LIB_OBJECT = $(BUILD_PATH)/librabbitbamtools.so $(BUILD_PATH)/librabbitbamread.so $(BUILD_PATH)/librabbitbamwrite.so

TARGET = $(BUILD_PATH)/rabbitbam

INCLUDE = /home/user_home/zz/tools/htslib-1.15.1/

LIB = /home/user_home/zz/tools/htslib-1.15.1/build/lib

SHARE = -fPIC -shared

api: block_mul.cpp
	$(CXX) $(CPPFLAGS) -I$(INCLUDE) -L$(LIB)     -I/home/user_home/zz/API/BamUniversalStatus -L/home/user_home/zz/API/BamUniversalStatus/build -lrabbitbamtools -lrabbitbamread -lrabbitbamwrite block_mul.cpp -o api

.PHONY:clean

clean:
	rm -rf api
