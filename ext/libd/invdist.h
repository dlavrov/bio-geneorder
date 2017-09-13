#ifndef INVDIST_H
#define INVDIST_H

#include "structs.h"

int invdist_noncircular ( Genome * g1, Genome * g2, int offset );
                          
int invdist_circular ( Genome * g1, Genome * g2);

int calculate_offset ( Genome * g1, Genome * g2);

void connected_component ( int size, distmem_t * distmem,
                           int * num_components );
                           
int num_cycles ( Genome * g1, Genome * g2 );

int num_breakpoints ( Genome * g1, Genome * g2 );

#endif
