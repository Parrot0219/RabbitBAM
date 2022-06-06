cmake_minimum_required(VERSION 3.17)
project(BamStatus)

set(CMAKE_CXX_STANDARD 14)

set(INC_DIR  /Users/zhaozhan/CLionProjects/BamUniversalStatus/htslib)
set(LINK_DIR /Users/zhaozhan/CLionProjects/BamUniversalStatus/htslib/build/lib)
include_directories(${INC_DIR})
link_directories(${LINK_DIR})

aux_source_directory(/Users/zhaozhan/CLionProjects/BamUniversalStatus/htslib HTSLIB)

add_executable(Config ${SAMTOOLS} ${HTSLIB} ${HTSLIB_TEXT}
        test/test_header.cpp main.cpp
        Duplicate.cpp Duplicate.h
        BamBlock.cpp BamBlock.h
        block_mul.cpp
        Buffer.cpp Buffer.h
        BamStatus.cpp BamStatus.h
        config.h
        Overrepresent.cpp Overrepresent.h
        test/test_sam.cpp test/test_sam_read.cpp
        blockingconcurrentqueue.h lightweightsemaphore.h concurrentqueue.h
        BamRead.h BamRead.cpp
        BamCompress.h BamCompress.cpp
        BamCompleteBlock.cpp BamCompleteBlock.h
        BamTools.cpp BamTools.h
        BamWrite.cpp BamWrite.h)


add_executable(BamStatus main.cpp)
add_executable(test_header test/test_header.cpp)
target_link_libraries(BamStatus libhts.a)