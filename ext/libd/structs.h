#ifndef STRUCTS_H
#define STRUCTS_H

#define ERR_CONTENT		-1
#define ERR_DUPLICATES	-2
#define ERR_CIRCULAR	-3
#define ERR_MULTICHR	-4

#define ERR_NOTIMPL		-5

#define HURDLE          1
#define GREATHURDLE (1<<1)
#define SUPERHURDLE (1<<2)

#include <stdlib.h>
#include <stdio.h>

typedef short intArray; /* T_ARRAY */

typedef struct {
	intArray * pi;
	bool circular;
	int len;
} Genome;

typedef struct component_struct
{
    int index;                  /* index of component's root */
    int oriented;               /* Boolean: Is component oriented */
    int blocks;                 /* Number of blocks in nonoriented component */
    int hurdle;                 /* Bitmask of HURDLE properties.
                                   Bit 0 = hurdle
                                   Bit 1 = wrap-around hurdle
                                   Bit 2 = superhurdle */
    int left;                   /* Index of component to the left of my rightmost block */
    int right;                  /* Index of component to the right of my rightmost block */
} component_t;

typedef struct distmem_struct
{
	int *perm1;                 /* DIST_INV: 2*num_genes + 2 */
    int *perm2;
    int *perm;
    int *done;
    int *greyEdges;

    int *stack;                 /* DIST_INV: size */
    int *oriented;
    int *cc;
    int *labeled;
    component_t *components;
    
    distmem_struct(int n){
    	perm1 = new int[n];
		perm2 = new int[n];
		perm = new int[n];
		done = new int[n];
		greyEdges = new int[n];
		stack = new int[n];
		oriented = new int[n];
		cc = new int[n];
		labeled = new int[n];
		components = new component_t[n];
	}
	~distmem_struct(){
		delete perm1;
		delete perm2;
		delete perm;
		delete done;
		delete greyEdges;
		delete stack;
		delete oriented;
		delete cc;
		delete labeled;
		delete components;
	}

} distmem_t;

#endif
