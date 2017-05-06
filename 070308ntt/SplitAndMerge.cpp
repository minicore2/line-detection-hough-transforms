#include "StdAfx.h"
#include "SplitAndMerge.h"

SplitAndMerge::SplitAndMerge(void)
{
}

SplitAndMerge::~SplitAndMerge(void)
{
}

void SplitAndMerge::SpaSMPInit(SpaSMPParams params) {  
	SPASMPPARAMS = params;    
	SpaCluInit(SPACLUPARAMS); 
}

void SplitAndMerge::SpaSMPCleanUp(void) {  
	SpaCluCleanUp(); 
}

int SplitAndMerge::SpaSMPExSeg2D(const SpaRawScan2D *scan, const int MaxNumOfSegs, int *NumOfSegs, SpaSegment2D *segments) {
  int i;
  SpaCluScan2D cluscan;
  SpaSegScan2D rawSegs;  

  SpaSegScan2DEmpty(scan->_cyclic, &rawSegs); // Initialization
  SpaCluExClusters2D(scan, &cluscan); // Forming clusters.

  for (i=0; i<cluscan._NClusters; ++i) // Extracting raw segments from the clusters.
    SpaSMPExSegs(&(cluscan._scan), cluscan._startIdx[i], cluscan._startIdx[i] + cluscan._length[i] - 1, 1, &rawSegs);

  if (SPASMPPARAMS.USEMERGEALGO) SpaMegMergeSeqSegs(&rawSegs); // Applying merging algorithm  
  SpaDelShortSegments(&rawSegs); // Removing short segments (always the last one)
  SpaSegsFromSegScan(&rawSegs, MaxNumOfSegs, NumOfSegs, segments); // Extracting results to the output 'segments'

  return SpaIntMin(MaxNumOfSegs, *NumOfSegs);
}
