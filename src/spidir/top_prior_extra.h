#ifndef SPIDIR_TOP_PRIOR_EXTRA_H
#define SPIDIR_TOP_PRIOR_EXTRA_H

#include "Tree.h"


namespace spidir {

extern "C" {

int inumHistories(int ngenes);

double numHistories(int ngenes);

int inumTopologyHistories(Tree *tree);

double numTopologyHistories(Tree *tree);


}

} // namespace spidir

#endif // SPIDIR_TOP_PRIOR_EXTRA_H
