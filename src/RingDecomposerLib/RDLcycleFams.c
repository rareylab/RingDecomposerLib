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

#include "RDLcycleFams.h"

#include <stdlib.h>
#include <limits.h>
#include <stdio.h>

#include "RDLapsp.h"
#include "RDLgraph.h"
#include "RDLutility.h"

/** returns 1 if the intersection of P(r,y) and P(r,z) is equal to {r};
0 otherwise */
static int RDL_pathsShareOnlyStart(unsigned r, unsigned y, unsigned z, RDL_graph *gra,
    RDL_sPathInfo *spi, unsigned *vertInRY, unsigned *vertInRZ)
{
  unsigned result = 0, i, pnt, count=0;
  /* vertices in P(r,y) and P(r,z); clear reusable scratch buffers */
  for(i=0; i<gra->V; ++i)
  {
    vertInRY[i] = 0;
    vertInRZ[i] = 0;
  }
  vertInRY[y] = 1;
  vertInRZ[z] = 1;
  /*find out which vertices are on the path P(r,y)*/
  pnt = y;
  do
  {
    pnt = spi->pred[r][pnt];
    vertInRY[pnt] = 1;
  }while(pnt != r);
  /*find out which vertices are on the path P(r,z)*/
  pnt = z;
  do
  {
    pnt = spi->pred[r][pnt];
    vertInRZ[pnt] = 1;
  }while(pnt != r);
  /*find out if more than one vertex is shared by the paths*/
  for(i=0; i<gra->V && count<2; ++i)
  {
    if(vertInRY[i] == 1 && vertInRZ[i] == 1)
    {
      ++count;
    }
  }
  /*if r is the only shared vertex of the two paths return 1*/
  if(count == 1 && (vertInRY[r] == 1) && (vertInRZ[r] == 1))
  {
    result = 1;
  }

  return result;
}

/** returns a cycle vector (element of {0,1}^m). odd (x=UINT_MAX) or even cycle.*/
static char *RDL_findPrototype(unsigned r, unsigned y, unsigned z, unsigned x,
    RDL_graph *gra, RDL_sPathInfo *spi)
{
  unsigned i, vert1, vert2;
  char *proto;

  proto = RDL_alloc(gra->E, sizeof(*proto));
  if (!proto) {
    return NULL;
  }
  for(i=0; i<gra->E; ++i)
  {
    proto[i] = 0;
  }
  /*path from r to y*/
  vert1 = y;
  do
  {
    vert2 = vert1;
    vert1 = spi->pred[r][vert1];
    proto[RDL_edgeId(gra, vert1, vert2)] = 1;
  }while(vert1 != r);
  /*path from r to z*/
  vert1 = z;
  do
  {
    vert2 = vert1;
    vert1 = spi->pred[r][vert1];
    proto[RDL_edgeId(gra, vert1, vert2)] = 1;
  }while(vert1 != r);
  if(x == UINT_MAX)/*odd cycle*/
  {
    proto[RDL_edgeId(gra,y,z)] = 1;
  }
  else /*even cycle*/
  {
    proto[RDL_edgeId(gra,y,x)] = 1;
    proto[RDL_edgeId(gra,z,x)] = 1;
  }
  return proto;
}

/** Adds an odd family when x == UINT_MAX, otherwise an even family. */
static int RDL_addFamily(unsigned r, unsigned y, unsigned z, unsigned x,
    RDL_graph *gra, RDL_sPathInfo *spi, RDL_cfURF *rc)
{
  RDL_cfam *family, **fams;
  if (rc->alloced == rc->nofFams) {
    if (rc->alloced > UINT_MAX / 2) {
      return 0;
    }
    fams = RDL_resize(rc->fams, (size_t)rc->alloced * 2, sizeof(*fams));
    if (!fams) {
      return 0;
    }
    rc->fams = fams;
    rc->alloced *= 2;
  }
  family = RDL_alloc(1, sizeof(*family));
  if (!family) {
    return 0;
  }
  family->prototype = RDL_findPrototype(r, y, z, x, gra, spi);
  if (!family->prototype) { RDL_free(family); return 0; }
  family->r = r; family->p = y; family->q = z; family->x = x;
  family->weight = spi->dist[r][y] + spi->dist[r][z] + (x == UINT_MAX ? 1 : 2);
  rc->fams[rc->nofFams++] = family;
  return 1;
}

/** finds a number of cycle families that contain at least all RELEVANT cycle
families - just like in Vismara's pseudocode */
static int RDL_vismara(RDL_cfURF *rc, RDL_graph *gra, RDL_sPathInfo *spi)
{
  unsigned i,j;
  int shared;
  unsigned rv,yv,zv,pv,qv; /*variables as in Vismara's algorithm, extended by a 'v'*/
  unsigned *evenCand; /*'S' in Vismara's algorithm*/
  unsigned *vertInRY = NULL, *vertInRZ = NULL;
  unsigned nofCandidates = 0;
  evenCand = RDL_alloc(gra->V, sizeof(*evenCand));
  vertInRY = RDL_alloc(gra->V, sizeof(*vertInRY));
  vertInRZ = RDL_alloc(gra->V, sizeof(*vertInRZ));
  if (!evenCand || !vertInRY || !vertInRZ) {
    RDL_free(evenCand);
    RDL_free(vertInRY);
    RDL_free(vertInRZ);
    return 0;
  }
  for(rv = 0; rv < gra->V; ++rv) {
    for(yv = 0; yv < gra->V; ++yv) {
      /*all yv reachable from rv respecting the ordering*/
      if(spi->reachable[rv][yv] == 1) {
        nofCandidates = 0;
        for(i = 0; i < gra->degree[yv]; ++i) {
          zv = gra->adjList[yv][i][0];
          /*all zv reachable from rv respecting the ordering and adjacent to yv*/
          if(spi->reachable[rv][zv] == 1) {
            if(spi->dist[rv][zv] + 1 == spi->dist[rv][yv]) {
              evenCand[nofCandidates] = zv;
              ++nofCandidates;
            }
            else if(spi->dist[rv][zv] != spi->dist[rv][yv] + 1
                && (gra->degree[zv] < gra->degree[yv] ||
                    (gra->degree[zv] == gra->degree[yv] && zv<yv))
                ) {
              /*add odd cycle rv-yv rv-zv zv-yv*/
              shared = RDL_pathsShareOnlyStart(rv, yv, zv, gra, spi,
                  vertInRY, vertInRZ);
              if (shared && !RDL_addFamily(rv, yv, zv, UINT_MAX, gra, spi, rc)) {
                goto fail;
              }
            }
            /*to fill dPaths in sPathInfo with the edges to r*/
            if(spi->dist[rv][zv] == 1) {
              if (!RDL_addEdge(spi->dPaths[rv], zv, rv)) {
                goto fail;
              }
            }
          }
        }
        /*any pair in evenCand*/
        for(i = 0; i < nofCandidates; ++i) {
          pv = evenCand[i];
          for(j = i+1; j < nofCandidates; ++j) {
            qv = evenCand[j];
            shared = RDL_pathsShareOnlyStart(rv, pv, qv, gra, spi,
                vertInRY, vertInRZ);
            if (shared) {
              /*add even cycle rv-pv rv-qv pv-yv-qv*/
              if (!RDL_addFamily(rv, pv, qv, yv, gra, spi, rc)) {
                goto fail;
              }
            }
          }
        }
        /*to fill dPaths in sPathInfo/fill U_r (see Vismara)*/
        for(i = 0; i < nofCandidates; ++i) {
          pv = evenCand[i];
          if (!RDL_addEdge(spi->dPaths[rv], yv, pv)) {
            goto fail;
          }
        }
      }
    }
  }

  RDL_free(evenCand);
  RDL_free(vertInRY);
  RDL_free(vertInRZ);
  return 1;
fail:
  RDL_free(evenCand);
  RDL_free(vertInRY);
  RDL_free(vertInRZ);
  return 0;
}

static int RDL_cycleFamsComp(const void *cf1, const void *cf2)
{
    if((*((RDL_cfam **)cf1))->weight < (*((RDL_cfam **)cf2))->weight)
    {
        return -1;
    }
    else if((*((RDL_cfam **)cf1))->weight > (*((RDL_cfam **)cf2))->weight)
    {
        return 1;
    } else {
      return 0;
    }
}

RDL_cfURF *RDL_findCycleFams(RDL_graph *gra, RDL_sPathInfo *spi)
{
  RDL_cfURF *rc = RDL_alloc(1, sizeof(*rc));
  if (!rc) {
    return NULL;
  }
  rc->nofFams = 0;
  rc->alloced = 64;
  rc->fams = RDL_alloc(rc->alloced, sizeof(*rc->fams));
  if (!rc->fams || !RDL_vismara(rc, gra, spi)) {
    RDL_deleteCycleFams(rc);
    return NULL;
  }

  /* sort by size */
  qsort(rc->fams, rc->nofFams, sizeof(*rc->fams), RDL_cycleFamsComp);
  return rc;
}

void RDL_deleteCycleFams(RDL_cfURF *rc)
{
  unsigned i;
  if (!rc) {
    return;
  }
  for(i=0; i<rc->nofFams; ++i)
  {
    RDL_free(rc->fams[i]->prototype);
    RDL_free(rc->fams[i]);
  }
  RDL_free(rc->fams);
  RDL_free(rc);
}

