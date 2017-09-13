#ifndef DISTANCES_H
#define DISTANCES_H

#include <math.h>
#include <vector>
#include <algorithm>
#include "structs.h"

std::vector<intArray> _adjacencies(Genome * pi, Genome * id);

int _breakpoints(Genome * pi, Genome * id);

int _inversions(Genome * pi, Genome * id);

int _DCJ(Genome * pi, Genome * id);

bool duplicates(Genome * pi);

bool unequal_content(Genome * pi, Genome * id);

#endif
