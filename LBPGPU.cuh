
#ifndef __LBPGPU_H__
#define __LBPGPU_H__

#include <iostream>

using namespace std;

typedef struct lbp_mapping {
	unsigned int * table;
	unsigned int samples;
	unsigned int number;
	float radius;
} LBPMapping;

void calcLBPGPU(const unsigned char * src, unsigned char * dst, const int width, const int height, const LBPMapping * map );


#endif
