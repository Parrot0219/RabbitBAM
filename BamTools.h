//
// Created by 赵展 on 2022/1/21.
//

#ifndef BAMSTATUS_BAMTOOLS_H
#define BAMSTATUS_BAMTOOLS_H

#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>

#include <htslib/khash.h>
#include "header.h"
#include <stdint.h>
#include "config.h"
#include <cstdlib>
#include <mutex>
#include <libdeflate.h>

#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8
//#define BGZF_MAX_BLOCK_SIZE 0x10000
#define BGZF_MAX_BLOCK_COMPLETE_SIZE 0x20000
typedef struct {
    int size;
    uint8_t *block;
    int64_t end_offset;
} cache_t;
KHASH_MAP_INIT_INT64(cache, cache_t)

typedef struct bam_block bam_block;
typedef struct bam_complete_block bam_complete_block;
struct bam_block{
    unsigned int errcode;
    unsigned char data[BGZF_MAX_BLOCK_SIZE];//0x1000
    unsigned int length;
    unsigned int pos;
    int64_t block_address;
    unsigned int split_pos;
    unsigned int bam_number;
};

struct bam_complete_block{
    unsigned int errcode;
    unsigned char *data = nullptr; // 指针变量未初始化会导致奇怪错误
    unsigned int data_size;
    unsigned int length;
    unsigned int pos;
    int64_t block_address;
};

struct bgzf_cache_t {
    khash_t(cache) *h;
    khint_t last_pos;
};
int sam_realloc_bam_data(bam1_t *b, size_t desired);
inline int realloc_bam_data(bam1_t *b, size_t desired);
inline int possibly_expand_bam_data(bam1_t *b, size_t bytes);
int bam_tag2cigar(bam1_t *b, int recal_bin, int give_warning); // return 0 if CIGAR is untouched; 1 if CIGAR is updated with CG
void bam_cigar2rqlens(int n_cigar, const uint32_t *cigar,
                      hts_pos_t *rlen, hts_pos_t *qlen);
void swap_data(const bam1_core_t *c, int l_data, uint8_t *data, int is_host);
inline int unpackInt16(const uint8_t *buffer);
int load_block_from_cache(BGZF *fp, int64_t block_address);
int check_header(const uint8_t *header);
const char *bgzf_zerr(int errnum, z_stream *zs);
int bgzf_uncompress(uint8_t *dst, size_t *dlen,const uint8_t *src, size_t slen,uint32_t expected_crc);
//int bgzf_uncompress(uint8_t *dst, unsigned int *dlen,const uint8_t *src, unsigned int slen,uint32_t expected_crc);
int fixup_missing_qname_nul(bam1_t *b) ;
int block_decode_func(struct bam_block *comp,struct bam_block *un_comp);
int read_block(BGZF *fp, struct bam_block *j);
int Rabbit_bgzf_read(struct bam_block *fq,void *data,unsigned int length);
void Rabbit_memcpy(void *target,unsigned char *source,unsigned int length);
int Rabbit_bgzf_read(struct bam_complete_block *fq,void *data,unsigned int length);
int read_bam(struct bam_block *fq,bam1_t *b,int is_be);
int read_bam(struct bam_complete_block *fq,bam1_t *b,int is_be);

/*
 *  获取bam——block块中最后一次完整分割的位置
 */
int find_divide_pos(bam_block *block,int last_pos=0);
int find_divide_pos(bam_complete_block *block,int last_pos=0);

/*
 *  获取bam——block块中最后一次完整分割的位置及其中存在的完整的bam数量
 *  返回值 ： 第一个为pos ，第二个为bam的数量
 */
std::pair<int,int> find_divide_pos_and_get_read_number(bam_block *block,int last_pos=0);
std::pair<int,int> find_divide_pos_and_get_read_number(bam_complete_block *block,int last_pos=0);

int change_data_size(bam_complete_block *block);


#endif //BAMSTATUS_BAMTOOLS_H
