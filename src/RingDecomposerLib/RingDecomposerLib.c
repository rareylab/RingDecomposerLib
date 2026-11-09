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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <string.h>
#include <float.h>
#include <math.h>

#include "RingDecomposerLib.h"

#include "RDLapsp.h"
#include "RDLbitset.h"
#include "RDLcycleFams.h"
#include "RDLdataStruct.h"
#include "RDLgraph.h"
#include "RDLhandler.h"
#include "RDLinfo.h"
#include "RDLrelation.h"
#include "RDLtarjan.h"
#include "RDLutility.h"

const unsigned RDL_RESERVED_START = 64;

/* define some constants */
const unsigned RDL_INVALID_RESULT = UINT_MAX;
const unsigned RDL_DUPLICATE_EDGE = UINT_MAX - 1;
const unsigned RDL_NO_RINGSYSTEM = UINT_MAX - 2;
const double RDL_INVALID_RC_COUNT = DBL_MAX;

/* public, set the global output function */
void RDL_setOutputFunction(RDL_outputFunction func)
{
  RDL_outputFunc = func ? func : RDL_writeNothing;
}

/* public */
RDL_graph *RDL_initNewGraph(unsigned V)
{
  return RDL_initNewGraph_g(V, 1);
}

/* public, add an edge with range check! */
unsigned RDL_addUEdge(RDL_graph *gra, RDL_node from, RDL_node to)
{
  return RDL_addUEdge_g(gra, from, to);
}

/* public, the calculation routine */
RDL_data *RDL_calculate(RDL_graph *gra)
{
  RDL_data *data;
  unsigned i, j, k, urf_index, rcf_index;
  unsigned nof_relevant_fams, nof_relevant_fams_sum;

  if (!gra) {
    RDL_outputFunc(RDL_ERROR, "The graph is NULL.\n");
    return NULL;
  }

  /* we can't calculate anything if there isn't at least one node */
  if(!gra->V) {
    RDL_outputFunc(RDL_ERROR, "The graph has no nodes.\n");
    return NULL;
  }

  data = RDL_alloc(1, sizeof(*data));

  if (!data) {
    return NULL;
  }

  /* FIRST STEP: TARJAN */
  data->bccGraphs = RDL_tarjanBCC(gra);
  if (!data->bccGraphs) {
    goto fail;
  }
  data->nofURFs = 0;
  data->nofRCFs = 0;

  /* allocate result structures */
  data->spiPerBCC = RDL_alloc(data->bccGraphs->nof_bcc, sizeof(*data->spiPerBCC));
  if (!data->spiPerBCC) {
    goto fail;
  }
  data->CFsPerBCC = RDL_alloc(data->bccGraphs->nof_bcc, sizeof(*data->CFsPerBCC));
  if (!data->CFsPerBCC) {
    goto fail;
  }
  data->urfInfoPerBCC = RDL_alloc(data->bccGraphs->nof_bcc, sizeof(*data->urfInfoPerBCC));
  if (!data->urfInfoPerBCC) {
    goto fail;
  }
  data->nofURFsPerBCC = RDL_alloc(data->bccGraphs->nof_bcc, sizeof(*data->nofURFsPerBCC));
  if (!data->nofURFsPerBCC) {
    goto fail;
  }
  data->nofRCFsPerBCC = RDL_alloc(data->bccGraphs->nof_bcc, sizeof(*data->nofRCFsPerBCC));
  if (!data->nofRCFsPerBCC) {
    goto fail;
  }

  /* the followeing steps are performed for each BCC */
  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    /* solve APSP problem */
    data->spiPerBCC[i] = RDL_AllPairsShortestPaths(data->bccGraphs->bcc_graphs[i]);
    if (!data->spiPerBCC[i]) {
      goto fail;
    }
    /* calculate RCFs with Vismara's algortihm */
    data->CFsPerBCC[i] = RDL_findCycleFams(data->bccGraphs->bcc_graphs[i], data->spiPerBCC[i]);
    if (!data->CFsPerBCC[i]) {
      goto fail;
    }
    if(data->CFsPerBCC[i]->nofFams > 0) {
      /* if there is at least one RCF, check URF relation */
      data->urfInfoPerBCC[i] = RDL_checkURFRelation(data->CFsPerBCC[i],
          data->bccGraphs->bcc_graphs[i], data->spiPerBCC[i]);
      if (!data->urfInfoPerBCC[i]) {
        goto fail;
      }
      data->nofURFsPerBCC[i] = data->urfInfoPerBCC[i]->nofURFs;

      nof_relevant_fams = 0;
      /* count RCFs */
      for (j = 0; j < data->CFsPerBCC[i]->nofFams; ++j) {
        if (data->CFsPerBCC[i]->fams[j]->mark) {
          ++nof_relevant_fams;
        }
      }

      nof_relevant_fams_sum = 0;
      for (j = 0; j < data->nofURFsPerBCC[i]; ++j) {
        nof_relevant_fams_sum += data->urfInfoPerBCC[i]->nofCFsPerURF[j];
      }

      if (nof_relevant_fams != nof_relevant_fams_sum) {
        RDL_outputFunc(RDL_ERROR, "different number of relevant families!\n");
        /* internal check, should never happen */
        assert(0);
      }
      data->nofRCFsPerBCC[i] = nof_relevant_fams;
    }
    else {
      data->nofURFsPerBCC[i] = 0;
      data->nofRCFsPerBCC[i] = 0;
    }
    data->nofURFs += data->nofURFsPerBCC[i];
    data->nofRCFs += data->nofRCFsPerBCC[i];
  }

  /* create a mapping from URFs to BCCs */
  data->urf_to_bcc = RDL_alloc(data->nofURFs, sizeof(*data->urf_to_bcc));
  if (!data->urf_to_bcc) {
    goto fail;
  }
  /* create a mapping from RCFs to URFs */
  data->rcf_to_urf = RDL_alloc(data->nofRCFs, sizeof(*data->rcf_to_urf));
  if (!data->rcf_to_urf) {
    goto fail;
  }
  urf_index = 0;
  rcf_index = 0;
  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    for (j = 0; j < data->nofURFsPerBCC[i]; ++j, ++urf_index) {
      data->urf_to_bcc[urf_index][0] = i;
      data->urf_to_bcc[urf_index][1] = j;
      for (k = 0; k < data->urfInfoPerBCC[i]->nofCFsPerURF[j]; ++k, ++rcf_index) {
        data->rcf_to_urf[rcf_index][0] = urf_index;
        data->rcf_to_urf[rcf_index][1] = k;
      }
    }
  }

  data->graph = gra;

  return data;
fail:
  /* Input ownership transfers only on success. */
  RDL_deleteData(data);
  return NULL;
}

/* public, free memory */
void RDL_deleteData(RDL_data *data)
{
  unsigned i;

  if (!data) {
    return;
  }
  for (i = 0; data->bccGraphs && i < data->bccGraphs->nof_bcc; ++i) {
    if (data->spiPerBCC) {
      RDL_deleteAPSP(data->spiPerBCC[i], data->bccGraphs->bcc_graphs[i]->V);
    }
    if (data->CFsPerBCC) {
      RDL_deleteCycleFams(data->CFsPerBCC[i]);
    }
    if (data->urfInfoPerBCC) {
      RDL_deleteURFInfo(data->urfInfoPerBCC[i]);
    }
  }
  RDL_free(data->spiPerBCC);
  RDL_free(data->CFsPerBCC);
  RDL_free(data->nofURFsPerBCC);
  RDL_free(data->nofRCFsPerBCC);
  RDL_free(data->urfInfoPerBCC);
  RDL_deleteBCCGraph(data->bccGraphs);
  RDL_free(data->urf_to_bcc);
  RDL_free(data->rcf_to_urf);

  RDL_deleteGraph(data->graph);
  RDL_free(data);
}

/* public, get the total number of URFs */
unsigned RDL_getNofURF(const RDL_data *data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }
  return data->nofURFs;
}

/* public, get the total number of RCFs */
unsigned RDL_getNofRCF(const RDL_data *data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }
  return data->nofRCFs;
}

/*
 * private
 * return the weight for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static unsigned RDL_getWeight_internal(const RDL_data *data,
    unsigned bcc_index, unsigned urf_internal_index, unsigned rcf_internal_index)
{
  return data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index][rcf_internal_index]->weight;
}

/* public, get weight for the given URF */
unsigned RDL_getWeightForURF(const RDL_data *data, unsigned index)
{
  unsigned bcc_index, internal_index;
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RESULT;
  }
  /* ger internal index */
  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];
  return RDL_getWeight_internal(data, bcc_index, internal_index, 0);
}

/* public, get weight for the given RCF */
unsigned RDL_getWeightForRCF(const RDL_data *data, unsigned index)
{
  unsigned bcc_index, urf_internal_index, rcf_internal_index, urf_index;
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RESULT;
  }

  /* map rcf to urf and urf to bcc */
  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];
  return RDL_getWeight_internal(data, bcc_index,
      urf_internal_index, rcf_internal_index);
}

/*
 * private
 * return the nodes for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static void RDL_getNodes_internal(const RDL_data *data, unsigned bcc_index,
    unsigned urf_internal_index, unsigned rcf_internal_index, char* atoms,
    char *visited)
{
  const RDL_cfam **URF;
  const RDL_graph* graph;

  graph = data->bccGraphs->bcc_graphs[bcc_index];

  URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];

  /* BFS on both shortest path graphs */
  RDL_giveVertices(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
      atoms, data->spiPerBCC[bcc_index], visited);

  memset(visited, 0, graph->V *  sizeof(*visited));
  RDL_giveVertices(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
      atoms, data->spiPerBCC[bcc_index], visited);

  if(URF[rcf_internal_index]->x < UINT_MAX) /*even cycle*/ {
    atoms[URF[rcf_internal_index]->x] = 1;
  }
}


/* private
 * gives an array of indices of nodes that are contained in the URF with the
 * given index. Array is terminated by UINT_MAX
 */
static RDL_node *RDL_getNodesURF(const RDL_data *data, unsigned index,
    unsigned *count)
{
  unsigned i, nofFams, bcc_index, internal_index;
  char *atoms = NULL, *visited = NULL;
  RDL_node *result;
  const RDL_graph* graph;

  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  atoms = RDL_alloc(graph->V, sizeof(*atoms));
  visited = RDL_alloc(graph->V, sizeof(*visited));
  if (!atoms || !visited) {
    goto fail;
  }

  nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
  /* iterate over all RCFs in this URFs and collect data in one bitset (atoms) */
  for(i=0; i<nofFams; ++i) {
    RDL_getNodes_internal(data, bcc_index, internal_index, i, atoms, visited);
  }

  result = RDL_bitset_to_array(atoms, graph->V,
      data->bccGraphs->node_from_bcc_mapping[bcc_index], count);
  RDL_free(atoms);
  RDL_free(visited);
  return result;
fail:
  RDL_free(atoms);
  RDL_free(visited);
  return NULL;
}

/* private
 * gives an array of indices of nodes that are contained in the RRF with the
 * given index. Array is terminated by UINT_MAX
 */
static RDL_node *RDL_getNodesRCF(const RDL_data *data, unsigned index,
    unsigned *count)
{
  unsigned bcc_index,
      urf_internal_index, rcf_internal_index, urf_index;
  char *atoms = NULL, *visited = NULL;
  RDL_node *result;
  const RDL_graph* graph;

  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  atoms = RDL_alloc(graph->V, sizeof(*atoms));
  visited = RDL_alloc(graph->V, sizeof(*visited));
  if (!atoms || !visited) {
    goto fail;
  }

  /* save the nodes into bitset */
  RDL_getNodes_internal(data, bcc_index, urf_internal_index, rcf_internal_index,
      atoms, visited);

  result = RDL_bitset_to_array(atoms, graph->V,
      data->bccGraphs->node_from_bcc_mapping[bcc_index], count);
  RDL_free(atoms);
  RDL_free(visited);
  return result;
fail:
  RDL_free(atoms);
  RDL_free(visited);
  return NULL;
}

/* public, calculate the number of RCs */
double RDL_getNofRC(const RDL_data *data)
{
  double result = 0, intermediate;
  unsigned i;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RC_COUNT;
  }

  /* add up results for URFs */
  for (i = 0; i < data->nofURFs; ++i) {
    intermediate = RDL_getNofRCForURF(data, i);
    if (intermediate == RDL_INVALID_RC_COUNT) {
      return RDL_INVALID_RC_COUNT;
    }

    result += intermediate;
  }

  return result;
}

/*
 * private
 * return the # RCs for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static double RDL_getNofRCForRCF_internal(const RDL_data *data,
    unsigned bcc_index, unsigned urf_internal_index,
    unsigned rcf_internal_index)
{
  double nofPaths1, nofPaths2, result;
  const double prod_limit = sqrt(DBL_MAX) - 1.0;
  const RDL_cfam** URF;
  const RDL_graph* graph;

  result = 0.0;

  graph = data->bccGraphs->bcc_graphs[bcc_index];

  URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];

  /* count the number of paths one both sides */
  nofPaths1 = RDL_countPaths(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
      graph->V, data->spiPerBCC[bcc_index]);
  nofPaths2 = RDL_countPaths(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
      graph->V, data->spiPerBCC[bcc_index]);

  if (nofPaths1 == RDL_INVALID_RC_COUNT || nofPaths2 == RDL_INVALID_RC_COUNT) {
    return RDL_INVALID_RC_COUNT;
  }
  result = nofPaths1 * nofPaths2;

  /* check if either of the number paths is larger than what we can multiply */
  if (nofPaths1 >= prod_limit || nofPaths2 >= prod_limit) {
    RDL_outputFunc(RDL_WARNING, "result overflow when counting paths!\n");
    return RDL_INVALID_RC_COUNT;
  }

  return result;
}

/* public, return the # RCs for given RCF */
double RDL_getNofRCForRCF(const RDL_data *data, unsigned index)
{
  unsigned bcc_index, urf_id, urf_internal_index, rcf_internal_index;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RC_COUNT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RC_COUNT;
  }

  /* map indices */
  urf_id = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];

  bcc_index = data->urf_to_bcc[urf_id][0];
  urf_internal_index = data->urf_to_bcc[urf_id][1];

  return RDL_getNofRCForRCF_internal(data, bcc_index,
      urf_internal_index, rcf_internal_index);
}

/* public, return the # RCs for given URF */
double RDL_getNofRCForURF(const RDL_data *data, unsigned index)
{
  unsigned i,nofFams, bcc_index, internal_index;
  double result=0, prod;

  const double sum_limit = DBL_MAX/2.0 - 1.0;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RC_COUNT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return RDL_INVALID_RC_COUNT;
  }

  /* map indices */
  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];

  nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];

  /* iterate all RCFs in this URF */
  for(i=0; i<nofFams; ++i) {
    prod = RDL_getNofRCForRCF_internal(data, bcc_index, internal_index, i);

    /* check if any summand is larger than limit */
    if (prod >= sum_limit || result >= sum_limit) {
      RDL_outputFunc(RDL_WARNING, "result overflow when counting paths!\n");
      return RDL_INVALID_RC_COUNT;
    }

    result += prod;
  }

  return result;
}

/* private
 * return the edges for INTERNAL indices,
 * urf_internal_index relative to bcc_index and rcf_internal_index relative to URF
 */
static void RDL_getEdges_internal(const RDL_data *data,
    unsigned bcc_index, unsigned urf_internal_index,
    unsigned rcf_internal_index, char* edges, char *visited)
{
  const RDL_cfam **URF;
  const RDL_graph* graph;

  graph = data->bccGraphs->bcc_graphs[bcc_index];

  URF = (const RDL_cfam **)data->urfInfoPerBCC[bcc_index]->URFs[urf_internal_index];

  /* perform BFS on both paths and collect edges */
  RDL_giveEdges(URF[rcf_internal_index]->r, URF[rcf_internal_index]->q,
      edges, graph, data->spiPerBCC[bcc_index], visited);

  memset(visited, 0, graph->V *  sizeof(*visited));
  RDL_giveEdges(URF[rcf_internal_index]->r, URF[rcf_internal_index]->p,
      edges, graph, data->spiPerBCC[bcc_index], visited);

  /*
   * in contrast to the nodes, we have to add additional edges, because we're
   * looking at cycles...
   */
  if(URF[rcf_internal_index]->x < UINT_MAX) /*even cycle*/ {
    edges[RDL_edgeId(graph,URF[rcf_internal_index]->q,URF[rcf_internal_index]->x)] = 1;
    edges[RDL_edgeId(graph,URF[rcf_internal_index]->p,URF[rcf_internal_index]->x)] = 1;
  }
  else /*odd cycle*/ {
    edges[RDL_edgeId(graph,URF[rcf_internal_index]->q,URF[rcf_internal_index]->p)] = 1;
  }
}


/* private
 * gives an array of indices of edges that are contained in the URF with the
 * given index. Array is terminated by UINT_MAX
 */
static unsigned *RDL_getEdgesURF(const RDL_data *data, unsigned index,
    unsigned *count)
{
  unsigned i, nofFams, bcc_index, internal_index;
  char *edges = NULL, *visited = NULL;
  unsigned *result;
  const RDL_graph* graph;

  /* map indices */
  bcc_index = data->urf_to_bcc[index][0];
  internal_index = data->urf_to_bcc[index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  edges = RDL_alloc(graph->E, sizeof(*edges));
  visited = RDL_alloc(graph->V, sizeof(*visited));
  if (!edges || !visited) {
    goto fail;
  }

  nofFams = data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];
  /* iterate over all RCFs in this URF and collect results in bitset (edges) */
  for(i=0; i<nofFams; ++i) {
    RDL_getEdges_internal(data, bcc_index, internal_index, i, edges, visited);
  }

  result = RDL_bitset_to_array(edges, graph->E,
      data->bccGraphs->edge_from_bcc_mapping[bcc_index], count);
  RDL_free(edges);
  RDL_free(visited);
  return result;
fail:
  RDL_free(edges);
  RDL_free(visited);
  return NULL;
}

/* private
 * gives an array of indices of edges that are contained in the RCF with the
 * given index. Array is terminated by UINT_MAX
 */
static unsigned *RDL_getEdgesRCF(const RDL_data *data, unsigned index,
    unsigned *count)
{
  unsigned urf_index, bcc_index,
      urf_internal_index, rcf_internal_index;
  char *edges = NULL, *visited = NULL;
  unsigned *result;
  const RDL_graph* graph;

  /* index mapping */
  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];
  graph = data->bccGraphs->bcc_graphs[bcc_index];
  edges = RDL_alloc(graph->E, sizeof(*edges));
  visited = RDL_alloc(graph->V, sizeof(*visited));
  if (!edges || !visited) {
    goto fail;
  }

  RDL_getEdges_internal(data, bcc_index, urf_internal_index,
      rcf_internal_index, edges, visited);

  result = RDL_bitset_to_array(edges, graph->E,
      data->bccGraphs->edge_from_bcc_mapping[bcc_index], count);
  RDL_free(edges);
  RDL_free(visited);
  return result;
fail:
  RDL_free(edges);
  RDL_free(visited);
  return NULL;
}


/* private, calls getNodesURF() or getEdgesURF() depending on mode 'a' or 'b' */
static unsigned *RDL_giveURF(const RDL_data *data, unsigned index, char mode)
{
  unsigned *result;
  if(mode == 'a')
  {
    result = RDL_getNodesURF(data, index, NULL);
  }
  else if(mode == 'b')
  {
    result = RDL_getEdgesURF(data, index, NULL);
  }
  else
  {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_giveURF()' with invalid mode '%c'\n", mode);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  return result;
}

/* private, calls getNodesRCF() or getEdgesRCF() depending on mode 'a' or 'b' */
static unsigned *RDL_giveRCF(const RDL_data *data, unsigned index, char mode)
{
  unsigned *result;
  if(mode == 'a')
  {
    result = RDL_getNodesRCF(data, index, NULL);
  }
  else if(mode == 'b')
  {
    result = RDL_getEdgesRCF(data, index, NULL);
  }
  else
  {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_giveRCF()' with invalid mode '%c'\n", mode);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  return result;
}

/* public, function for retrieving nodes in an URF */
unsigned RDL_getNodesForURF(const RDL_data *data, unsigned index, RDL_node **ptr)
{
  unsigned count;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  (*ptr) = RDL_getNodesURF(data, index, &count);

  if (!*ptr) {
    return RDL_INVALID_RESULT;
  }
  return count;
}

/* public, function for retrieving nodes in an RCF */
unsigned RDL_getNodesForRCF(const RDL_data *data, unsigned index, RDL_node **ptr)
{
  unsigned count;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  (*ptr) = RDL_getNodesRCF(data, index, &count);

  if (!*ptr) {
    return RDL_INVALID_RESULT;
  }
  return count;
}

static RDL_edge *RDL_mapEdgeIndices(const RDL_data *data,
    const unsigned *edgeIndices, unsigned count)
{
  unsigned i;
  RDL_edge *result;

  result = RDL_alloc(count, sizeof(*result));
  if (!result) {
    return NULL;
  }
  for (i = 0; i < count; ++i) {
    result[i][0] = data->graph->edges[edgeIndices[i]][0];
    result[i][1] = data->graph->edges[edgeIndices[i]][1];
  }
  return result;
}

/* public, function for retrieving edges of URF */
unsigned RDL_getEdgesForURF(const RDL_data *data, unsigned index, RDL_edge **ptr)
{
  unsigned count;
  unsigned *edgeIndices = NULL;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  edgeIndices = RDL_getEdgesURF(data, index, &count);
  if (!edgeIndices) {
    goto fail;
  }
  *ptr = RDL_mapEdgeIndices(data, edgeIndices, count);
  RDL_free(edgeIndices);
  edgeIndices = NULL;
  if (!*ptr) {
    goto fail;
  }
  return count;
fail:
  RDL_free(edgeIndices);
  *ptr = NULL;
  return RDL_INVALID_RESULT;
}

/* public, function for retrieving edges of RCF */
unsigned RDL_getEdgesForRCF(const RDL_data *data, unsigned index, RDL_edge **ptr)
{
  unsigned count;
  unsigned *edgeIndices = NULL;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  edgeIndices = RDL_getEdgesRCF(data, index, &count);
  if (!edgeIndices) {
    goto fail;
  }
  *ptr = RDL_mapEdgeIndices(data, edgeIndices, count);
  RDL_free(edgeIndices);
  edgeIndices = NULL;
  if (!*ptr) {
    goto fail;
  }
  return count;
fail:
  RDL_free(edgeIndices);
  *ptr = NULL;
  return RDL_INVALID_RESULT;
}

/* public, to delete result of RDL_getCyclesChar() */
void RDL_deleteCyclesChar(char **cycles)
{
  unsigned i;
  for(i=0; cycles[i]!=NULL; ++i)
  {
    RDL_free(cycles[i]);
  }
  RDL_free(cycles);
}

/* public, delete cycle */
void RDL_deleteCycle(RDL_cycle *cycle)
{
  if (!cycle) {
    return;
  }
  RDL_free(cycle->edges);
  RDL_free(cycle);
}

/* public, delete cycles */
void RDL_deleteCycles(RDL_cycle **cycles, unsigned number)
{
  unsigned i;
  if (!cycles) {
    return;
  }
  for(i=0; i<number; ++i)
  {
    RDL_deleteCycle(cycles[i]);
  }
  RDL_free(cycles);
}

/*
 * private
 * Gives all families containing the object, which can be an atom or a bond.
mode:
  - 'a': nodes
  - 'b': edges
returns an array of integers containing all indices of families containing the
object. The array ends with a terminating UINT_MAX on its last position and has
to be deallocated with 'RDL_free()'
*/
static unsigned *RDL_listFamilies(const RDL_data *data, unsigned object,
    char mode, char family, unsigned *count)
{
  unsigned *result = NULL;
  void *resized;
  unsigned nof_families;
  unsigned nextfree=0, alloced=RDL_RESERVED_START;
  unsigned *objects;
  unsigned i,j, bcc_index, target_bcc, urf_index;
  char contained, can_be_contained;

  if(!(mode == 'a' || mode == 'b'))
  {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_listFamilies()' with invalid mode '%c'\n",mode);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  unsigned* (*giveFamily) (const RDL_data*, unsigned, char) = NULL;

  if (family == 'u') {
    nof_families = data->nofURFs;
    giveFamily = RDL_giveURF;
  }
  else if (family == 'r') {
    nof_families = data->nofRCFs;
    giveFamily = RDL_giveRCF;
  }
  else {
    RDL_outputFunc(RDL_ERROR, "tried to call 'RDL_listFamilies()' with invalid family '%c'\n",family);
    /* cannot occur when using interface */
    assert(0);
    return NULL;
  }

  result = RDL_alloc((size_t)alloced + 1, sizeof(*result));

  if (!result) {
    goto fail;
  }
  for(i = 0; i < nof_families; ++i)
  {
    if (family == 'u') {
      urf_index = i;
    }
    else if (family == 'r') {
      urf_index = data->rcf_to_urf[i][0];
    }

    bcc_index = data->urf_to_bcc[urf_index][0];

    /*
     * check if object belongs to this BCC at all,
     * if not it can't be in the URF either
     */
    can_be_contained = 0;
    if (mode == 'a') {
      for (j = 0; j < data->bccGraphs->nof_bcc_per_node[object]; ++j) {
        target_bcc = data->bccGraphs->node_to_bcc_mapping[object][2*j];
        if (target_bcc == bcc_index) {
          can_be_contained = 1;
          break;
        }
      }
    }
    else {
      target_bcc = data->bccGraphs->edge_to_bcc_mapping[object][0];
      if (target_bcc == bcc_index) {
        can_be_contained = 1;
      }
    }

    contained = 0;
    if (can_be_contained) {
      objects = giveFamily(data, i, mode);
      if (!objects) {
        goto fail;
      }
      for(j=0; objects[j]<UINT_MAX; ++j)
      {
        if(objects[j] == object)
        {
          contained = 1;
          break;
        }
      }
      RDL_free(objects);
    }
    if (contained) {
      if(nextfree == alloced) {
        if (alloced > UINT_MAX / 2) {
          goto fail;
        }
        alloced *= 2;
        resized = RDL_resize(result, (size_t)alloced + 1, sizeof(*result));
        if (!resized) {
          goto fail;
        }
        result = resized;
      }
      result[nextfree++] = i;
    }
  }

  result[nextfree] = UINT_MAX;
  *count = nextfree;
  return result;
fail:
  RDL_free(result);
  return NULL;
}

/* private, get families containing a node */
static unsigned RDL_getFamiliesContainingNode_internal(
    const RDL_data *data, RDL_node object, unsigned **ptr,
    char family)
{
  unsigned count;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (object >= data->graph->V) {
    RDL_outputFunc(RDL_ERROR, "invalid node: %u\n", object);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = NULL;
    return 0;
  }

  *ptr = RDL_listFamilies(data, object, 'a', family, &count);

  if (!*ptr) {
    return RDL_INVALID_RESULT;
  }
  return count;
}

/* public, get the RCF indices with this node */
unsigned RDL_getRCFsContainingNode(const RDL_data *data, RDL_node object,
                                   unsigned **ptr)
{
  return RDL_getFamiliesContainingNode_internal(data, object, ptr, 'r');
}

/* public, get the URF indices with this node */
unsigned RDL_getURFsContainingNode(const RDL_data *data, RDL_node object,
                                   unsigned **ptr)
{
  return RDL_getFamiliesContainingNode_internal(data, object, ptr, 'u');
}

/* private, get the family indices with this edge */
static unsigned RDL_getFamiliesContainingEdge_internal(
    const RDL_data *data, RDL_node node1, RDL_node node2, unsigned **ptr,
    char family)
{
  unsigned count, edge;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (node1 >= data->graph->V || node2 >= data->graph->V) {
    RDL_outputFunc(RDL_ERROR, "invalid edge: %u %u\n", node1, node2);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  edge = RDL_edgeId(data->graph, node1, node2);

  /* check if this edge exists at all */
  if (edge == RDL_INVALID_RESULT) {
    RDL_outputFunc(RDL_ERROR, "invalid edge: %u %u\n", node1, node2);
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = NULL;
    return 0;
  }

  *ptr = RDL_listFamilies(data, edge, 'b', family, &count);

  if (!*ptr) {
    return RDL_INVALID_RESULT;
  }
  return count;
}


/* public, get the RCF indices with this edge */
unsigned RDL_getRCFsContainingEdge(const RDL_data *data, RDL_node node1,
                                   RDL_node node2, unsigned **ptr)
{
  return RDL_getFamiliesContainingEdge_internal(data, node1, node2, ptr, 'r');
}

/* public, get the URF indices with this edge */
unsigned RDL_getURFsContainingEdge(const RDL_data *data, RDL_node node1,
                                   RDL_node node2, unsigned **ptr)
{
  return RDL_getFamiliesContainingEdge_internal(data, node1, node2, ptr, 'u');
}

/* public, get the number of URFs that this node is part of */
unsigned RDL_getNofURFContainingNode(const RDL_data *data, RDL_node node)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getURFsContainingNode(data, node, &arr);
  RDL_free(arr);
  return number;
}

/* public, get the number of URFs that this edge is part of */
unsigned RDL_getNofURFContainingEdge(const RDL_data *data, RDL_node node1, RDL_node node2)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getURFsContainingEdge(data, node1, node2, &arr);
  RDL_free(arr);
  return number;
}

/* public, get the number of RCFs that this node is part of */
unsigned RDL_getNofRCFContainingNode(const RDL_data *data, RDL_node node)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getRCFsContainingNode(data, node, &arr);
  RDL_free(arr);
  return number;
}

/* public, get the number of RCFs that this edge is part of */
unsigned RDL_getNofRCFContainingEdge(const RDL_data *data, RDL_node node1, RDL_node node2)
{
  unsigned number;
  unsigned *arr;
  number = RDL_getRCFsContainingEdge(data, node1, node2, &arr);
  RDL_free(arr);
  return number;
}

/* construct a SSSR */
unsigned RDL_getSSSR(const RDL_data *data, RDL_cycle ***ptr)
{
  unsigned j, k, rcf=0, urf=0;
  RDL_cycle **result = NULL;
  void *resized;
  unsigned size, added, currBond;
  unsigned bcc_index = 0, internal_index;
  char* prototype;
  const RDL_graph* curr_graph;
  unsigned total_counter = 0, alloced;
  unsigned success = 0;
  unsigned edge;
  unsigned char *currentCycle;
  unsigned char **basisCycles = NULL;
  unsigned char **prototypesForGaussian = NULL;
  unsigned currentBasisCyclesSize = 0;
  unsigned oldBasisCyclesSize;
  unsigned compressedSize = 0;
  const unsigned char* empty_cycle = NULL;


  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = NULL;
    return 0;
  }

  total_counter = 0;
  if (data->graph->E < data->graph->V) {
    /*
     * if E - V + 1 < 1  <=> E - V < 0 <=> E < V
     *
     * graph is not connected and our heuristic for reservation
     * fails => allocate something ;)
     * */
    alloced = RDL_RESERVED_START;
  }
  else {
    alloced = data->graph->E - data->graph->V + 1;
  }
  result = RDL_alloc(alloced, sizeof(*result));
  if (!result) {
    goto fail;
  }

  /* calculate basis for each BCC (and therefore for each CC) */

  for (bcc_index = 0; bcc_index < data->bccGraphs->nof_bcc; ++bcc_index) {
    curr_graph = data->bccGraphs->bcc_graphs[bcc_index];

    size=curr_graph->E - curr_graph->V + 1;

    basisCycles = RDL_alloc(size, sizeof(*basisCycles));

    if (!basisCycles) {
      goto fail;
    }
    prototypesForGaussian = RDL_alloc(data->nofURFsPerBCC[bcc_index], sizeof(*prototypesForGaussian));
    if (!prototypesForGaussian) {
      goto fail;
    }

    /* copy the prototypes for gaussian elimination */
    /*
     * it is sufficient to take the first prototype for each URF
     * (the other's are linearly dependent anyway
     */
    for (internal_index = 0;
        internal_index < data->nofURFsPerBCC[bcc_index]; ++internal_index) {
      compressedSize = RDL_bitset_compressed(&(prototypesForGaussian[internal_index]),
          data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->prototype,
          curr_graph->E);
      if (compressedSize == RDL_INVALID_RESULT) {
        goto fail;
      }
    }

    empty_cycle = RDL_alloc(compressedSize, sizeof(*empty_cycle));

    if (!empty_cycle) {
      goto fail;
    }

    added = 0;
    currentBasisCyclesSize = 0;

    for (internal_index = 0;
        internal_index < data->nofURFsPerBCC[bcc_index]; ++internal_index, ++urf) {

      currentCycle = RDL_alloc(compressedSize, sizeof(*currentCycle));

      if (!currentCycle) {
        goto fail;
      }
      memcpy(currentCycle, prototypesForGaussian[internal_index],
          compressedSize * sizeof(*currentCycle));

      /*
       * basisCycles is in row echelon form!
       * [this corresponds to "add row operation"]
       */
      for (k = 0; k < currentBasisCyclesSize; ++k) {
        if (RDL_bitset_test(currentCycle, k)) {
          RDL_bitset_xor_inplace(currentCycle, basisCycles[k], compressedSize);
        }
      }

      if (RDL_bitset_empty(currentCycle, empty_cycle, compressedSize)) {
        RDL_free(currentCycle);
        continue;
      }
      /*
       * the current cycle is indepedent of equal and smaller
       * sized cycles => part of the basis
       */
      oldBasisCyclesSize = currentBasisCyclesSize;
      basisCycles[currentBasisCyclesSize] = currentCycle;
      ++currentBasisCyclesSize;

      /*
       * this is where the magic is happening: ensure that the basis is
       * indeed in row echelon form (we use the property throughout the algorithm)
       *
       * swap columns such that indeed the current ring with index currentBasisCycleSize-1
       * has a 1 at column currentBasisCycleSize-1
       * => ROW ECHELON FORM
       *
       * [this corresponds to the "swap columns operation" of gaussian elimination]
       */
      if (!RDL_bitset_test(currentCycle, oldBasisCyclesSize)) {
        for (k = oldBasisCyclesSize + 1; k < curr_graph->E; ++k) {
          if (RDL_bitset_test(currentCycle, k)) {
            RDL_swap_columns(basisCycles, currentBasisCyclesSize, oldBasisCyclesSize, k);
            RDL_swap_columns(prototypesForGaussian, data->nofURFsPerBCC[bcc_index], oldBasisCyclesSize, k);
            break;
          }
        }
      }

      prototype = data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->prototype;

      /* can happen for unconnected graphs */
      if (total_counter >= alloced) {
        if (alloced > UINT_MAX / 2) {
          goto fail;
        }
        alloced *= 2;
        resized = RDL_resize(result, alloced, sizeof(*result));
        if (!resized) {
          goto fail;
        }
        result = resized;
      }

      currBond = 0;
      result[total_counter] = RDL_alloc(1, sizeof(**result));
      if (!result[total_counter]) {
        goto fail;
      }
      result[total_counter]->edges = RDL_alloc(data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->weight, sizeof(*(result[total_counter]->edges)));
      if (!result[total_counter]->edges) { RDL_deleteCycle(result[total_counter]); goto fail; }
      result[total_counter]->weight = data->urfInfoPerBCC[bcc_index]->URFs[internal_index][0]->weight;
      result[total_counter]->urf = urf;
      result[total_counter]->rcf = rcf;
      for(j=0; j < curr_graph->E; ++j) {
        if(prototype[j] == 1) {
          edge = data->bccGraphs->edge_from_bcc_mapping[bcc_index][j];
          result[total_counter]->edges[currBond][0] = data->graph->edges[edge][0];
          result[total_counter]->edges[currBond][1] = data->graph->edges[edge][1];
          ++currBond;
        }
      }

      rcf += data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[internal_index];

      /*if enough inRDL_dependent cycles were found, break out of the loop*/
      ++total_counter;
      if(++added == size) {
        break;
      }
    }

    for (k = 0; k < currentBasisCyclesSize; ++k) {
      RDL_free(basisCycles[k]);
    }
    RDL_free(basisCycles); basisCycles = NULL; currentBasisCyclesSize = 0;

    for (k = 0; k < data->nofURFsPerBCC[bcc_index]; ++k) {
      RDL_free(prototypesForGaussian[k]);
    }
    RDL_free(prototypesForGaussian); prototypesForGaussian = NULL;
    RDL_free((void*)empty_cycle); empty_cycle = NULL;
  }



  success = 1;
fail:
  if (basisCycles) {
    for (k = 0; k < currentBasisCyclesSize; ++k) {
      RDL_free(basisCycles[k]);
    }
  }
  RDL_free(basisCycles);
  if (prototypesForGaussian) {
    for (k = 0; k < data->nofURFsPerBCC[bcc_index]; ++k) {
      RDL_free(prototypesForGaussian[k]);
    }
  }
  RDL_free(prototypesForGaussian); RDL_free((void*)empty_cycle);
  if (!success) {
    RDL_deleteCycles(result, total_counter);
    *ptr = NULL;
    return RDL_INVALID_RESULT;
  }
  *ptr = result;
  return total_counter;
}

/* public, get the RCPs */
unsigned RDL_getRCPrototypes(const RDL_data *data, RDL_cycle ***ptr)
{
  RDL_cycle **result = NULL;
  unsigned nof_urf;
  unsigned nofRel=0, rcf=0;
  unsigned currFam=0, currEdge;
  unsigned i,j,k;
  const char* prototype;
  unsigned bcc_index, mapped_edge;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (data->nofURFs < 1) {
    (*ptr) = NULL;
    return 0;
  }

  for (i = 0; i < data->bccGraphs->nof_bcc; ++i) {
    for (j = 0; j < data->nofURFsPerBCC[i]; ++j) {
      nofRel += data->urfInfoPerBCC[i]->nofCFsPerURF[j];
    }
  }

  /*allocate space*/
  result = RDL_alloc(nofRel, sizeof(*result));
  if (!result) {
    goto fail;
  }

  for (bcc_index = 0; bcc_index < data->bccGraphs->nof_bcc; ++bcc_index) {
    nof_urf = data->nofURFsPerBCC[bcc_index];
    for (i = 0; i < nof_urf; ++i) {
      for (j = 0; j < data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[i]; ++j, ++rcf) {
        prototype = data->urfInfoPerBCC[bcc_index]->URFs[i][j]->prototype;
        result[currFam] = RDL_alloc(1, sizeof(**result));
        if (!result[currFam]) {
          goto fail;
        }
        result[currFam]->edges = RDL_alloc(data->urfInfoPerBCC[bcc_index]->URFs[i][j]->weight, sizeof(*result[currFam]->edges));
        if (!result[currFam]->edges) { RDL_deleteCycle(result[currFam]); goto fail; }
        result[currFam]->weight = data->urfInfoPerBCC[bcc_index]->URFs[i][j]->weight;
        result[currFam]->urf = i;
        result[currFam]->rcf = rcf;
        currEdge = 0;
        for (k = 0; k < data->bccGraphs->bcc_graphs[bcc_index]->E; ++k) {
          if (prototype[k] == 1) {
            mapped_edge = data->bccGraphs->edge_from_bcc_mapping[bcc_index][k];
            result[currFam]->edges[currEdge][0] = data->graph->edges[mapped_edge][0];
            result[currFam]->edges[currEdge][1] = data->graph->edges[mapped_edge][1];
            ++currEdge;
          }
        }
        ++currFam;
      }
    }
  }

  (*ptr) = result;
  return nofRel;
fail:
  RDL_deleteCycles(result, currFam);
  *ptr = NULL;
  return RDL_INVALID_RESULT;
}

/* public, get an RDL_cycleIterator for RCFs */
RDL_cycleIterator* RDL_getRCyclesForRCFIterator(const RDL_data *data, unsigned index)
{
  RDL_cycleIterator* it;
  unsigned rcf_internal_index, urf_internal_index, urf_index;
  unsigned bcc_index;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return NULL;
  }

  if (index >= data->nofRCFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return NULL;
  }

  urf_index = data->rcf_to_urf[index][0];
  rcf_internal_index = data->rcf_to_urf[index][1];
  bcc_index = data->urf_to_bcc[urf_index][0];
  urf_internal_index = data->urf_to_bcc[urf_index][1];

  it = RDL_initCycleIterator(RDL_RCF_IT,
      rcf_internal_index, rcf_internal_index,
      urf_internal_index, urf_internal_index,
      bcc_index, bcc_index,
      'b', data);

  return it;
}

/* public, get an RDL_cycleIterator for URFs */
RDL_cycleIterator* RDL_getRCyclesForURFIterator(const RDL_data *data, unsigned index)
{
  RDL_cycleIterator* it;
  unsigned urf_internal_index, bcc_index;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return NULL;
  }

  if (index >= data->nofURFs) {
    RDL_outputFunc(RDL_ERROR, "invalid index: %u\n", index);
    return NULL;
  }

  bcc_index = data->urf_to_bcc[index][0];
  urf_internal_index = data->urf_to_bcc[index][1];

  it = RDL_initCycleIterator(
      RDL_URF_IT,
      0, 0,
      urf_internal_index, urf_internal_index,
      bcc_index, bcc_index,
      'b', data);

  return it;
}

/* public, get an RDL_cycleIterator for all RCs */
RDL_cycleIterator* RDL_getRCyclesIterator(const RDL_data *data)
{
  RDL_cycleIterator* it;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return NULL;
  }

  it = RDL_initCycleIterator(
      RDL_ALL_IT,
      0, 0,
      0, 0,
      0, data->bccGraphs->nof_bcc-1,
      'b', data);

  return it;
}

/* private, get RCs from an RDL_cycleIterator */
static unsigned RDL_getRCycles_internal(
    const RDL_data *data, RDL_cycle ***ptr, RDL_cycleIterator* it)
{
  RDL_cycle **result = NULL;
  void *resized;
  unsigned alloced, nextfree;

  alloced = RDL_RESERVED_START;
  nextfree = 0;
  result = RDL_alloc(alloced, sizeof(*result));
  if (!result) {
    goto fail;
  }

  while (!RDL_cycleIteratorAtEnd(it)) {
    if(nextfree == alloced) {
      if (alloced > UINT_MAX / 2) {
        goto fail;
      }
      alloced *= 2;
      resized = RDL_resize(result, alloced, sizeof(*result));
      if (!resized) {
        goto fail;
      }
      result = resized;
    }
    result[nextfree] = RDL_cycleIteratorGetCycle(it);
    if (!result[nextfree]) {
      goto fail;
    }
    ++nextfree;

    if (!RDL_cycleIteratorNext(it)) {
      goto fail;
    }
  }
  /* Shrink the returned array to its exact size. */
  {
    RDL_cycle **shrunk = RDL_resize(result, nextfree, sizeof(*result));
    if (!shrunk) {
      goto fail;
    }
    result = shrunk;
  }

  RDL_deleteCycleIterator(it);
  *ptr = result;
  return nextfree;
fail:
  RDL_deleteCycles(result, nextfree);
  RDL_deleteCycleIterator(it);
  *ptr = NULL;
  return RDL_INVALID_RESULT;
}

/* public, get all RCs */
unsigned RDL_getRCycles(const RDL_data *data, RDL_cycle ***ptr)
{
  RDL_cycleIterator* it;

  it = RDL_getRCyclesIterator(data);

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  return RDL_getRCycles_internal(data, ptr, it);
}

/* public, get all RCs for given URF */
unsigned RDL_getRCyclesForURF(const RDL_data *data, unsigned index, RDL_cycle ***ptr)
{
  RDL_cycleIterator* it;

  it = RDL_getRCyclesForURFIterator(data, index);

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  return RDL_getRCycles_internal(data, ptr, it);
}

/* public, get all RCs for given RCF */
unsigned RDL_getRCyclesForRCF(const RDL_data *data, unsigned index, RDL_cycle ***ptr)
{
  RDL_cycleIterator* it;

  it = RDL_getRCyclesForRCFIterator(data, index);

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  return RDL_getRCycles_internal(data, ptr, it);
}

/* public, utils function for converting cycle to bitset representation */
unsigned RDL_translateCycArray(const RDL_data *data, RDL_cycle **array,
                               unsigned number, char ***ptr)
{
  unsigned i,j,RDL_edgeIdx;
  char **result;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  if (number < 1) {
    (*ptr) = NULL;
    return 0;
  }

  result = RDL_alloc2DCharArray(number, data->graph->E);
  if (!result) { *ptr = NULL; return RDL_INVALID_RESULT; }
  for(i=0; i<number;  ++i)/*initialize to 0s*/
  {
    for(j=0; j<data->graph->E; ++j)
    {
      result[i][j] = 0;
    }
  }

  for(i=0; i<number; ++i)
  {
    for(j=0; j<array[i]->weight; ++j)
    {
      RDL_edgeIdx = RDL_edgeId(data->graph, array[i]->edges[j][0],
                               array[i]->edges[j][1]);
      result[i][RDL_edgeIdx] = 1;
    }
  }
  (*ptr) = result;
  return number;
}

/* public, delete an edge array */
void RDL_deleteEdgeIdxArray(char **cycles, unsigned num)
{
  RDL_delete2DCharArray(cycles, num);
}

/* public, get the edges of the graph */
unsigned RDL_getEdgeArray(const RDL_data *data, RDL_edge **ptr)
{
  RDL_edge *result;
  unsigned idx;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*ptr) = NULL;
    return RDL_INVALID_RESULT;
  }

  result = RDL_alloc(data->graph->E, sizeof(*result));

  if (!result) {
    goto fail;
  }
  for(idx=0; idx<data->graph->E; ++idx)
  {
    result[idx][0] = data->graph->edges[idx][0];
    result[idx][1] = data->graph->edges[idx][1];
  }
  (*ptr) = result;
  return data->graph->E;

fail:
  *ptr = NULL;
  return RDL_INVALID_RESULT;
}

/* public, get the number of edges in the graph */
unsigned RDL_getNofEdges(const RDL_data *data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }
  return data->graph->E;
}

/* public, get the number of ringsystems calculated */
unsigned RDL_getNofRingsystems(const RDL_data* data)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->nof_bcc;
}

/* public, get the number of nodes for the ring sytem */
unsigned RDL_getNofNodesForRingsystem(const RDL_data *data, unsigned idx)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->bcc_graphs[idx]->V;
}

/* public, get the number of edges for the ring sytem */
unsigned RDL_getNofEdgesForRingsystem(const RDL_data *data, unsigned idx)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->bcc_graphs[idx]->E;

}

/* public, get the edge for the ring sytem */
unsigned RDL_getEdgesForRingsystem(
    const RDL_data *data, unsigned idx, RDL_edge** edges)
{
  unsigned i, edge;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*edges) = NULL;
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    (*edges) = NULL;
    return RDL_INVALID_RESULT;
  }

  (*edges) = RDL_alloc(data->bccGraphs->bcc_graphs[idx]->E, sizeof(**edges));

  if (!(*edges)) {
    goto fail;
  }

  for (i = 0; i < data->bccGraphs->bcc_graphs[idx]->E; ++i) {
    edge = data->bccGraphs->edge_from_bcc_mapping[idx][i];
    (*edges)[i][0] = data->graph->edges[edge][0];
    (*edges)[i][1] = data->graph->edges[edge][1];
  }

  return data->bccGraphs->bcc_graphs[idx]->E;

fail:
  *edges = NULL;
  return RDL_INVALID_RESULT;
}

/* public, get the nodes for the ring sytem */
unsigned RDL_getNodesForRingsystem(
    const RDL_data *data, unsigned idx, RDL_node** nodes)
{
  unsigned i, node;

  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    (*nodes) = NULL;
    return RDL_INVALID_RESULT;
  }

  if(idx >= data->bccGraphs->nof_bcc) {
    RDL_outputFunc(RDL_ERROR, "idx %d is out of range!\n", idx);
    (*nodes) = NULL;
    return RDL_INVALID_RESULT;
  }

  (*nodes) = RDL_alloc(data->bccGraphs->bcc_graphs[idx]->V, sizeof(**nodes));

  if (!(*nodes)) {
    goto fail;
  }

  for (i = 0; i < data->bccGraphs->bcc_graphs[idx]->V; ++i) {
    node = data->bccGraphs->node_from_bcc_mapping[idx][i];
    (*nodes)[i] = node;
  }

  return data->bccGraphs->bcc_graphs[idx]->V;

fail:
  *nodes = NULL;
  return RDL_INVALID_RESULT;
}

/* public, get the ringsystem an edge belongs to */
unsigned RDL_getRingsystemForEdge(
    const RDL_data* data, unsigned from, unsigned to)
{
  unsigned edge;

  edge = RDL_getEdgeId(data, from, to);

  if (edge == RDL_INVALID_RESULT) {
    return RDL_INVALID_RESULT;
  }

  return data->bccGraphs->edge_to_bcc_mapping[edge][0];
}

/* public, get the ID of an edge */
unsigned RDL_getEdgeId(const RDL_data *data, unsigned from, unsigned to)
{
  if (!data) {
    RDL_outputFunc(RDL_ERROR, "RDL_data is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (!data->graph) {
    RDL_outputFunc(RDL_ERROR, "RDL_graph is NULL!\n");
    return RDL_INVALID_RESULT;
  }

  if (from >= data->graph->V || to >= data->graph->V) {
    RDL_outputFunc(RDL_ERROR, "invalid edge %u %u\n", from, to);
    return RDL_INVALID_RESULT;
  }

  return RDL_edgeId(data->graph, from, to);
}
