/*
 * This file is part of the RingDecomposerLib, licensed
 * under BSD New license (see LICENSE in the root directory).
 * Copyright (c) 2016
 * University of Hamburg, ZBH - Center for Bioinformatics
 * Niek Andresen, Florian Flachsenberg, Matthias Rarey
 * 
 * Please cite:
 * 
 * Kolodzik, A.; Urbaczek, S.; Rarey, M.
 * Unique Ring Families: A Chemically Meaningful Description
 * of Molecular Ring Topologies.
 * J. Chem. Inf. Model., 2012, 52 (8), pp 2013-2021
 * 
 * Flachsenberg, F.; Andresen, N.; Rarey, M.
 * RingDecomposerLib: An Open-Source Implementation of
 * Unique Ring Families and Other Cycle Bases.
 * J. Chem. Inf. Model., 2017, 57 (2), pp 122-126
 */

#include "RDLinfo.h"

#include <stdio.h>
#include <stdlib.h>
#include "RDLcycleFams.h"
#include "RDLrelation.h"
#include "RDLutility.h"

/** allocates the matrices in urfInfo */
RDL_URFinfo *RDL_initUrfInfo(RDL_cfURF *RCFs)
{
  RDL_URFinfo *info = RDL_alloc(1, sizeof(*info));
  unsigned i, w = 0;
  if (!info) {
    return NULL;
  }
  info->nofWeights = 1;
  for (i = 1; i < RCFs->nofFams; ++i) {
    if (RCFs->fams[i]->weight != RCFs->fams[i - 1]->weight) {
      ++info->nofWeights;
    }
  }
  info->nofProtos = RDL_alloc(info->nofWeights, sizeof(*info->nofProtos));
  info->URFrel = RDL_alloc(info->nofWeights, sizeof(*info->URFrel));
  if (!info->nofProtos || !info->URFrel) {
    goto fail;
  }
  for (i = 0; i < RCFs->nofFams; ++i) {
    if (i && RCFs->fams[i]->weight != RCFs->fams[i - 1]->weight) {
      ++w;
    }
    ++info->nofProtos[w];
  }
  for (i = 0; i < info->nofWeights; ++i) {
    info->URFrel[i] = RDL_alloc2DCharArray(info->nofProtos[i], info->nofProtos[i]);
    if (!info->URFrel[i]) {
      goto fail;
    }
  }
  return info;
fail:
  RDL_deleteURFInfo(info);
  return NULL;
}

void RDL_deleteURFInfo(RDL_URFinfo *uInfo)
{
  unsigned i;
  if (!uInfo) {
    return;
  }
  for(i=0; uInfo->URFrel && i<uInfo->nofWeights; ++i)
  {
    RDL_delete2DCharArray(uInfo->URFrel[i], uInfo->nofProtos ? uInfo->nofProtos[i] : 0);
  }
  RDL_free(uInfo->URFrel);
  RDL_free(uInfo->nofProtos);
  for(i=0; uInfo->URFs && i<uInfo->nofURFs; ++i)
  {
    RDL_free(uInfo->URFs[i]);
  }
  RDL_free(uInfo->URFs);
  RDL_free(uInfo->nofCFsPerURF);
  RDL_free(uInfo);
}

static char **RDL_allocProtoFlags(RDL_URFinfo *info)
{
  unsigned i;
  char **flags = RDL_alloc(info->nofWeights, sizeof(*flags));
  if (!flags) {
    return NULL;
  }
  for (i = 0; i < info->nofWeights; ++i) {
    flags[i] = RDL_alloc(info->nofProtos[i], sizeof(**flags));
    if (!flags[i]) { RDL_delete2DCharArray(flags, info->nofWeights); return NULL; }
  }
  return flags;
}

/** returns the number of URFs. */
unsigned RDL_countURFs(RDL_URFinfo *uInfo)
{
  unsigned URFRelCount=0, weightIdx, i, j;
  char **alreadyInURF; /*for each weight: an array storing if a RCF is already
                        part of a URF*/
  char increase;

  alreadyInURF = RDL_allocProtoFlags(uInfo);
  if (!alreadyInURF) {
    return RDL_INVALID_RESULT;
  }
  for(i=0; i<uInfo->nofWeights; ++i)
  {
    for(j=0; j<uInfo->nofProtos[i]; ++j)
    {
      alreadyInURF[i][j] = 0;
    }
  }

  /*count number of 1s indicating URF-relation*/
  for(weightIdx=0; weightIdx<uInfo->nofWeights; ++weightIdx)
  {
    for(i=0; i<uInfo->nofProtos[weightIdx]; ++i)
    {
      if(alreadyInURF[weightIdx][i] == 1 || uInfo->URFrel[weightIdx][i][i] == 0)
      {
        continue;
      }
      increase = 0;
      for(j=i; j<uInfo->nofProtos[weightIdx]; ++j)
      {
        if(uInfo->URFrel[weightIdx][i][j] == 1)
        {
          alreadyInURF[weightIdx][j] = 1;
          if(increase == 0)
          {
            ++URFRelCount;
            increase = 1;
          }
        }
      }
    }
  }

  for(i=0; i<uInfo->nofWeights; ++i)
  {
    RDL_free(alreadyInURF[i]);
  }
  RDL_free(alreadyInURF);
  return URFRelCount;
}

/** returns the index in the array of RCFs (fams) that the RCF with the weight
at the index "weight" and position j has */
unsigned RDL_idxWeight(RDL_URFinfo *uInfo, unsigned weight, unsigned j)
{
  unsigned i, sum = 0;
  for(i=0; i<weight; ++i)
  {
    sum += uInfo->nofProtos[i];
  }
  return sum + j;
}

/** adds the CF with the index 'RCFIdx' to the URF with the index 'URFIdx'.
sNeeds the number of RCFs already in the URF as third parameter
'previouslyAlloced' */
int RDL_addRCFtoURF(unsigned RCFIdx, unsigned URFIdx, unsigned previouslyAlloced,
                     RDL_cfam ***URFs, RDL_cfURF *CFs)
{
  RDL_cfam **families = RDL_resize(URFs[URFIdx], (size_t)previouslyAlloced + 1,
                                  sizeof(*families));
  if (!families) {
    return 0;
  }
  URFs[URFIdx] = families;
  families[previouslyAlloced] = CFs->fams[RCFIdx];
  return 1;
}

int RDL_fillURFs(RDL_URFinfo *uInfo, RDL_cfURF *CFs)
{
  RDL_cfam ***URFs;
  unsigned *nofRCFs;
  unsigned currURFIdx=0;
  unsigned i,j,k;
  char **alreadyInURF; /*for each weight: an array storing if a RCF is already
                        part of a URF*/

  alreadyInURF = RDL_allocProtoFlags(uInfo);
  if (!alreadyInURF) {
    return 0;
  }
  for(i=0; i<uInfo->nofWeights; ++i)
  {
    for(j=0; j<uInfo->nofProtos[i]; ++j)
    {
      alreadyInURF[i][j] = 0;
    }
  }

  URFs = RDL_alloc(uInfo->nofURFs, sizeof(*URFs));
  nofRCFs = RDL_alloc(uInfo->nofURFs, sizeof(*nofRCFs));
  uInfo->URFs = URFs;
  uInfo->nofCFsPerURF = nofRCFs;
  if (!URFs || !nofRCFs) {
    goto fail;
  }
  for(i=0; i<uInfo->nofURFs; ++i)
  {
    nofRCFs[i] = 0;
  }

  for(i=0; i<uInfo->nofWeights; ++i)
  {
    for(j=0; j<uInfo->nofProtos[i]; ++j)
    {
      /*if this family is not already part of a URF*/
      if((uInfo->URFrel[i][j][j] == 1) && (alreadyInURF[i][j] == 0))
      {
        for(k=j; k<uInfo->nofProtos[i]; ++k)
        {
          if(uInfo->URFrel[i][j][k] == 1)
          {
            /*add to current URF the CFs with the given Idx*/
            if (!RDL_addRCFtoURF(RDL_idxWeight(uInfo, i, k), currURFIdx, nofRCFs[currURFIdx], URFs,
                                 CFs)) {
              goto fail;
            }
            ++nofRCFs[currURFIdx];
            alreadyInURF[i][k] = 1;
          }
        }
        ++currURFIdx;
      }
    }
  }

  for(i=0; i<uInfo->nofWeights; ++i)
  {
    RDL_free(alreadyInURF[i]);
  }
  RDL_free(alreadyInURF);
  return 1;
fail:
  RDL_delete2DCharArray(alreadyInURF, uInfo->nofWeights);
  return 0;
}

RDL_URFinfo *RDL_checkURFRelation(RDL_cfURF *RCFs, RDL_graph *graph, RDL_sPathInfo* spi)
{
  RDL_URFinfo *uInfo = RDL_initUrfInfo(RCFs);
  if (!uInfo) {
    return NULL;
  }
  if (!RDL_findRelations(RCFs, graph, uInfo, spi)) {
    goto fail;
  }

  uInfo->nofURFs = RDL_countURFs(uInfo);
  if (uInfo->nofURFs == RDL_INVALID_RESULT) { uInfo->nofURFs = 0; goto fail; }
  if (!RDL_fillURFs(uInfo, RCFs)) {
    goto fail;
  }
  return uInfo;
fail:
  RDL_deleteURFInfo(uInfo);
  return NULL;
}

