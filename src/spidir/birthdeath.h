/*=============================================================================

  Matt Rasmussen
  Copyright 2010-2011

  Birth-death models

=============================================================================*/


#ifndef SPIDIR_BIRTHDEATH_H
#define SPIDIR_BIRTHDEATH_H


namespace spidir {

extern "C" {


double birthDeathCount(int ngenes, float time, float birth, float death);
double birthDeathCounts(int start, int end, float time, 
                        float birth, float death);
double birthDeathCountsLog(int start, int end, float time, 
                           float birth, float death);
double sampleBirthWaitTime1(float T, float birth, float death);


}

} // namespace spidir

#endif // SPIDIR_BIRTHDEATH_H
