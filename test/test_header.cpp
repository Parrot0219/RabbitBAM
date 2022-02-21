//
// Created by 赵展 on 2021/6/11.
//
#include <htslib/sam.h>
#include <htslib/bgzf.h>
#include <htslib/hfile.h>
#include <zlib.h>
#include <htslib/khash.h>
#include "header.h"
#include "config.h"
#include <iostream>
#define BLOCK_HEADER_LENGTH 18
#define BLOCK_FOOTER_LENGTH 8
#define TDEF(x_)
#define TSTART(x_)
#define TEND(x_)
#define TPRINT(x_, str)
using namespace std;

int main(){
    TDEF(fq)
    TSTART(fq)
    samFile *sin;
    sam_hdr_t *hdr;
    bam1_t *b;
    if ((sin=sam_open("/home/old_home/haoz/workspace/data/NC/NC_T_1.sorted.bam", "r"))==NULL){
        printf("Can`t open this file!\n");
        return 0;
    }
    if ((hdr = sam_hdr_read(sin)) == NULL) {
        return  0;
    }
    // HD
    kstring_t *tmp =new kstring_t ;

    sam_hdr_find_tag_hd(hdr,"SO",tmp);
    cout << "SO : " << tmp->s << endl;
    sam_hdr_find_tag_id(hdr,"RG",NULL,NULL,"PL",tmp);
    cout << "PL ： " << tmp->s << endl;

    // CHR
    // sam_hdr_tid2name tid -> name
    for (int i=0;i<hdr->n_targets;i++){
        cout << "Tid -> chr : " << sam_hdr_tid2name(hdr,i) << endl;
    }

    cout << "chrX -> tid  :" << sam_hdr_name2tid(hdr,"chrX") << endl;
    cout << "chrY -> tid  :" << sam_hdr_name2tid(hdr,"chrY") << endl;

    sam_close(sin);
    TEND(fq)
    TPRINT(fq,"change time is : ");
}

