//
// Created by 赵展 on 2022/1/21.
//

#include "BamRead.h"

BamRead::BamRead(int BufferSize){
    readBlockSize = BufferSize+1;
    readBlock = new bam_block*[readBlockSize];
    read_bg = 0;read_ed = BufferSize-1;
    for (int i=read_bg;i<=read_ed;i++) readBlock[i] = new bam_block;

    consumerBlockSize=BufferSize+5;
    consumerBlock = new bam_block*[consumerBlockSize];
    consumer_bg = 1;consumer_ed = 0;
    blockNum = 0;

    read_complete = false;
}

bam_block* BamRead::getEmpty() {
    while ((read_ed+1)%readBlockSize == read_bg){
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
    }
    int num = read_bg;
    read_bg=(read_bg+1)%readBlockSize;
    return readBlock[num];
}

void BamRead::inputBlock(bam_block* block){
    consumerBlock[(consumer_ed+1)%consumerBlockSize] = block;
    consumer_ed = (consumer_ed+1)%consumerBlockSize;
}

std::pair<bam_block*, int> BamRead::getReadBlock(){
    mtx_consumer.lock();
    while ((consumer_ed+1)%consumerBlockSize == consumer_bg){
        mtx_consumer.unlock();
        std::this_thread::sleep_for(std::chrono::milliseconds(5));
        if (read_complete && (consumer_ed+1)%consumerBlockSize == consumer_bg ) return std::pair<bam_block *,int>(NULL,-1);
        mtx_consumer.lock();
    }
    int num_point = consumer_bg;
    consumer_bg = (consumer_bg+1)%consumerBlockSize;
    int num_block = blockNum++;
    mtx_consumer.unlock();
    return std::pair<bam_block*, int>(consumerBlock[num_point],num_block);
}

void BamRead::backBlock(bam_block* block){
    mtx_read.lock();
    readBlock[(read_ed+1)%readBlockSize] = block;
    read_ed = (read_ed+1)%readBlockSize;
    mtx_read.unlock();
}

void BamRead::ReadComplete(){
    read_complete=true;
}