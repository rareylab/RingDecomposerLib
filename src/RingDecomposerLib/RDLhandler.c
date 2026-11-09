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

#include "RDLhandler.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "RDLapsp.h"
#include "RDLbitset.h"
#include "RDLgraph.h"
#include "RDLstack.h"
#include "RDLcycleFams.h"
#include "RDLinfo.h"
#include "RDLdataStruct.h"
#include "RDLtarjan.h"
#include "RDLutility.h"

static void RDL_findIndegree(unsigned a, unsigned b,
    unsigned *indegree, const RDL_sPathInfo *spi)
{
/* similar to the function "List_Paths" from Vismara*/
  unsigned i, vertex;

  indegree[b] = 0;

  if(a==b)
  {
    return;
  }
  /*for each vertex adjacent to b in U_a*/
  for(i=0; i<spi->dPaths[a]->degree[b]; ++i)
  {
    vertex = spi->dPaths[a]->adjList[b][i][0];
    if (indegree[vertex] == UINT_MAX) {
      RDL_findIndegree(a, vertex, indegree, spi);
    }
    indegree[vertex] += 1;
  }
}


/* count the number of paths from a to b */
double RDL_countPaths(unsigned a, unsigned b,
    unsigned V, const RDL_sPathInfo *spi)
{
  unsigned i, curr, vertex, empty_end=0;
  unsigned *indegree = NULL, *empty_set = NULL;
  double *nof_paths = NULL, result;

  nof_paths = RDL_alloc(V, sizeof(*nof_paths));

  if (!nof_paths) {
    goto fail;
  }
  indegree = RDL_alloc(V, sizeof(*indegree));
  if (!indegree) {
    goto fail;
  }
  empty_set = RDL_alloc(V, sizeof(*empty_set));
  if (!empty_set) {
    goto fail;
  }
  for (i = 0; i < V; ++i) {
    nof_paths[i] = 0.0;
    indegree[i] = UINT_MAX;
  }

  RDL_findIndegree(a, b, indegree, spi);

  for (i = 0; i < V; ++i) {
    if (indegree[i] == 0) {
      empty_set[empty_end] = i;
      ++empty_end;
      nof_paths[i] = 1.0;
    }
  }

  if (empty_end != 1 || empty_set[0] != b) {
    RDL_outputFunc(RDL_ERROR, "invalid topological sort!");
    /* should never happen */
    assert(0);
    RDL_free(nof_paths);
    RDL_free(indegree);
    RDL_free(empty_set);
    return DBL_MAX;
  }

  /* look at the nodes in an topological sort */
  while (empty_end) {
    curr = empty_set[empty_end-1];
    --empty_end;
    for(i=0; i<spi->dPaths[a]->degree[curr]; ++i) {
      vertex = spi->dPaths[a]->adjList[curr][i][0];
      --indegree[vertex];
      nof_paths[vertex] += nof_paths[curr];
      if (!indegree[vertex]) {
        empty_set[empty_end] = vertex;
        ++empty_end;
      }
    }
  }

  result = nof_paths[a];

  RDL_free(nof_paths);
  RDL_free(indegree);
  RDL_free(empty_set);

  return result;
fail:
  RDL_free(nof_paths); RDL_free(indegree); RDL_free(empty_set);
  return RDL_INVALID_RC_COUNT;
}


void RDL_giveVertices(unsigned a, unsigned b, char *array,
    const RDL_sPathInfo *spi, char *visited)
{
 /* similar to the function "List_Paths" from Vismara*/
  unsigned i, vertex;

  visited[b] = 1;

  if(a==b)
  {
    array[a] = 1;
    return;
  }
  array[b] = 1;
  /*for each vertex adjacent to b in U_a*/
  for(i=0; i<spi->dPaths[a]->degree[b]; ++i)
  {
    vertex = spi->dPaths[a]->adjList[b][i][0];
    if (!visited[vertex]) {
      RDL_giveVertices(a, vertex, array, spi, visited);
    }
  }
}

void RDL_giveEdges(unsigned a, unsigned b, char *array,
    const RDL_graph *gra, const RDL_sPathInfo *spi, char* visited)
{
/* similar to the function "List_Paths" from Vismara*/
  unsigned i, vertex, edge;

  if(a==b) {
    return;
  }

  visited[b] = 1;

  /*for each vertex adjacent to b in U_a*/
  for(i=0; i<spi->dPaths[a]->degree[b]; ++i)
  {
    vertex = spi->dPaths[a]->adjList[b][i][0];
    edge = RDL_edgeId(gra,b,vertex);
    array[edge] = 1;
    if (!visited[vertex]) {
      RDL_giveEdges(a, vertex, array, gra, spi, visited);
    }
  }
}

/* this struct contains the COMPLETE state for path enumeration */
typedef struct RDL_pathIterator {
  RDL_stack *stack;
  unsigned char *compressed_path;
  unsigned compressed_number;
  char mode;
  const RDL_graph *gra;
  const RDL_sPathInfo *spi;
  char end;
} RDL_pathIterator;

typedef struct RDL_lpStackElement {
    unsigned r;
    unsigned p;
    unsigned ai_index;
    unsigned char *split_copy;
} RDL_lpStackElement;

static void RDL_lpStackElement_delete(RDL_lpStackElement* element)
{
  if (element->split_copy) {
    RDL_free(element->split_copy);
  }
  RDL_free(element);
}

static RDL_pathIterator* RDL_pathIteratorNext(RDL_pathIterator *it);
static void RDL_deletePathIterator(RDL_pathIterator *it);

/** just like the function "List_Paths" by Vismara: finds all shortest paths
between r and p recursively and writes them into the array "paths". 'mode'
determines if the paths are stored as arrays of vertices or of edges.*/
static RDL_pathIterator* RDL_listPaths(unsigned r, unsigned p,
                   char mode, const RDL_graph *gra, const RDL_sPathInfo *spi)
{
  unsigned number;
  RDL_lpStackElement *new_element = NULL;

  RDL_pathIterator *it = RDL_alloc(1, sizeof(*it));
  if (!it) {
    return NULL;
  }

  it->stack = RDL_stack_new();
  if (!it->stack) {
    goto fail;
  }
  it->mode = mode;

  if(mode == 'a') {
    number = gra->V;
  }
  else {
    number = gra->E;
  }
  it->compressed_number = RDL_bitset_init(&it->compressed_path, number);
  if (it->compressed_number == RDL_INVALID_RESULT) {
    goto fail;
  }

  it->spi = spi;
  it->gra = gra;
  it->end = 0;
  new_element = RDL_alloc(1, sizeof(*new_element));
  if (!new_element) {
    goto fail;
  }
  new_element->r = r;
  new_element->p = p;
  new_element->ai_index = 0;
  new_element->split_copy = NULL;

  if (!RDL_stack_push(it->stack, new_element)) {
    goto fail;
  }
  new_element = NULL;

  if (!RDL_pathIteratorNext(it)) {
    goto fail;
  }

  return it;
fail:
  RDL_free(new_element);
  RDL_deletePathIterator(it);
  return NULL;
}

/* free path iterator */
static void RDL_deletePathIterator(RDL_pathIterator *it)
{
  RDL_lpStackElement *top;
  if (!it) {
    return;
  }
  if (it->compressed_path) {
    RDL_free(it->compressed_path);
  }
  /* stack can be non-empty if premature terminating iteration */
  while (it->stack && !RDL_stack_empty(it->stack)) {
    top = RDL_stack_top(it->stack);
    RDL_lpStackElement_delete(top);
    RDL_stack_pop(it->stack);
  }
  RDL_stack_delete(it->stack);
  RDL_free(it);
}

/* move on to the next path */
static RDL_pathIterator* RDL_pathIteratorNext(RDL_pathIterator *it)
{
  unsigned vertex;
  RDL_lpStackElement *new_element = NULL, *top;

  if (it->end) {
    RDL_outputFunc(RDL_ERROR, "You tried to next an ended iterator!\n");
    return NULL;
  }

  if (RDL_stack_empty(it->stack)) {
    it->end = 1;
  }

  while (!RDL_stack_empty(it->stack)) {
    top = RDL_stack_top(it->stack);

    if(it->mode == 'a') {
      RDL_bitset_set(it->compressed_path, top->p);
    }

    if(top->r == top->p) {
      RDL_stack_pop(it->stack);
      RDL_lpStackElement_delete(top);
      break;
    }

    if(it->spi->dPaths[top->r]->degree[top->p] > 1)/*split*/
    {
      vertex = it->spi->dPaths[top->r]->adjList[top->p][top->ai_index][0];

      if (top->ai_index == 0) {
        top->split_copy = RDL_alloc(it->compressed_number, sizeof(*top->split_copy));
        if (!top->split_copy) {
          goto fail;
        }
        /*make a copy of the path so far*/
        memcpy(top->split_copy, it->compressed_path,
            it->compressed_number * sizeof(*top->split_copy));
      }
      else /*not continuing the last path*/
      {
        /*get the path up until the split*/
        memcpy(it->compressed_path, top->split_copy,
            it->compressed_number * sizeof(*top->split_copy));
      }
      if(it->mode == 'b')
      {
        RDL_bitset_set(it->compressed_path, RDL_edgeId(it->gra,top->p,vertex));
      }

      new_element = RDL_alloc(1, sizeof(*new_element));

      if (!new_element) {
        goto fail;
      }
      new_element->r = top->r;
      new_element->p = vertex;
      new_element->ai_index = 0;
      new_element->split_copy = NULL;
      ++top->ai_index;
      if (top->ai_index >= it->spi->dPaths[top->r]->degree[top->p]) {
        RDL_stack_pop(it->stack);
        RDL_lpStackElement_delete(top);
      }

      if (!RDL_stack_push(it->stack, new_element)) {
        goto fail;
      }
      new_element = NULL;
    }
    else/*not a split*/
    {
      vertex = it->spi->dPaths[top->r]->adjList[top->p][0][0];
      if(it->mode == 'b')
      {
        RDL_bitset_set(it->compressed_path, RDL_edgeId(it->gra,top->p,vertex));
      }
      top->p = vertex;
      top->ai_index = 0;
    }
  }

  return it;
fail:
  RDL_free(new_element);
  it->end = 1;
  return NULL;
}

/* this iterator containes the COMPLETE state for RC enumeration */
struct RDL_cycleIterator {
    RDL_pathIterator *it1;
    RDL_pathIterator *it2;
    char mode;
    char end;
    unsigned char *compressed_cycle;
    const RDL_data *data;
    unsigned rcf_index;
    unsigned rcf_index_max;
    unsigned urf_index;
    unsigned urf_index_max;
    unsigned bcc_index;
    unsigned bcc_index_max;
    unsigned running_rcf;
    unsigned running_urf;
    RDL_IteratorType type;
};

/* combine to paths to a relevant cycle */
static void RDL_combinePathsToCycle(unsigned char* compressed_cycle, const RDL_pathIterator *it1,
    const RDL_pathIterator *it2, const RDL_cfam* cfam)
{
  memcpy(compressed_cycle, it1->compressed_path,
      it1->compressed_number * sizeof(*compressed_cycle));
  RDL_bitset_or_inplace(compressed_cycle, it2->compressed_path,
      it2->compressed_number);
  if(it1->mode == 'a') {
    /*even cycle*/
    if(cfam->x < UINT_MAX) {
      RDL_bitset_set(compressed_cycle, cfam->x);
    }
  }
  else {
    if(cfam->x < UINT_MAX) {
      /*even cycle*/
      RDL_bitset_set(compressed_cycle, RDL_edgeId(it1->gra,cfam->p,cfam->x));
      RDL_bitset_set(compressed_cycle, RDL_edgeId(it1->gra,cfam->q,cfam->x));
    }
    else {
      /*odd cycle*/
      RDL_bitset_set(compressed_cycle, RDL_edgeId(it1->gra,cfam->p,cfam->q));
    }
  }
}

/* move on to the next relevant cycle */
RDL_cycleIterator* RDL_cycleIteratorNext(RDL_cycleIterator* it)
{
  const RDL_URFinfo *uInfo;
  const RDL_cfam *cfam;
  const RDL_graph *gra;
  const RDL_sPathInfo *spi;
  unsigned number;

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    return NULL;
  }

  if (RDL_cycleIteratorAtEnd(it)) {
    RDL_outputFunc(RDL_ERROR, "Cannot advance iterator at end!\n");
    return NULL;
  }

  uInfo = it->data->urfInfoPerBCC[it->bcc_index];
  cfam = uInfo->URFs[it->urf_index][it->rcf_index];
  gra = it->data->bccGraphs->bcc_graphs[it->bcc_index];
  spi = it->data->spiPerBCC[it->bcc_index];

  /* first call */
  if (it->it2 == NULL) {
    it->it2 = RDL_listPaths(cfam->r, cfam->q, it->mode, gra, spi);
    if (!it->it2) {
      goto fail;
    }
  }
  else {
    /* advance second iterator */
    if (!RDL_pathIteratorNext(it->it2)) {
      goto fail;
    }
  }

  if (it->it1 == NULL) {
    it->it1 = RDL_listPaths(cfam->r, cfam->p, it->mode, gra, spi);
    if (!it->it1) {
      goto fail;
    }
  }

  /* if second iterator is at end, increment first */
  if (it->it2->end) {
    /* do not increment first, if not possible! */
    if (!it->it1->end) {
      /* advance first iterator */
      if (!RDL_pathIteratorNext(it->it1)) {
        goto fail;
      }

      /* re initialize second iterator */
      RDL_deletePathIterator(it->it2);
      it->it2 = RDL_listPaths(cfam->r, cfam->q, it->mode, gra, spi);
      if (!it->it2) {
        goto fail;
      }
    }
  }

  /* if any of the iterators is STILL at the end, we have to change families */
  if (it->it1->end || it->it2->end) {
    /* if RCFs are at end */
    if (it->rcf_index == it->rcf_index_max) {
      if (it->urf_index == it->urf_index_max) {
        if (it->bcc_index == it->bcc_index_max) {
          it->end = 1;
        }
        else {
          ++it->bcc_index;
          it->urf_index = 0;
          it->urf_index_max = it->data->nofURFsPerBCC[it->bcc_index] - 1;
          it->rcf_index = 0;
          it->rcf_index_max = it->data->urfInfoPerBCC[it->bcc_index]->nofCFsPerURF[it->urf_index] - 1;
          if (it->mode == 'a') {
            number = it->data->bccGraphs->bcc_graphs[it->bcc_index]->V;
          }
          else {
            number = it->data->bccGraphs->bcc_graphs[it->bcc_index]->E;
          }
          RDL_free(it->compressed_cycle);
          it->compressed_cycle = NULL;
          if (RDL_bitset_init(&it->compressed_cycle, number) == RDL_INVALID_RESULT) {
            goto fail;
          }
          /* reset all iterators */
          RDL_deletePathIterator(it->it1);
          it->it1 = NULL;
          RDL_deletePathIterator(it->it2);
          it->it2 = NULL;
          /* now try again with new RCF */
          if (!RDL_cycleIteratorNext(it)) {
            goto fail;
          }
        }
      }
      else {
        ++it->urf_index;
        it->rcf_index = 0;
        it->rcf_index_max = uInfo->nofCFsPerURF[it->urf_index] - 1;
        /* reset all iterators */
        RDL_deletePathIterator(it->it1);
        it->it1 = NULL;
        RDL_deletePathIterator(it->it2);
        it->it2 = NULL;
        /* now try again with new RCF */
        if (!RDL_cycleIteratorNext(it)) {
          goto fail;
        }
      }
    }
    else {
      ++it->rcf_index;
      ++it->running_rcf;
      /* reset all iterators */
      RDL_deletePathIterator(it->it1);
      it->it1 = NULL;
      RDL_deletePathIterator(it->it2);
      it->it2 = NULL;
      /* now try again with new RCF */
      if (!RDL_cycleIteratorNext(it)) {
        goto fail;
      }
    }
  }
  else {
    /* otherwise, a new cycle is born ... */
    RDL_combinePathsToCycle(it->compressed_cycle, it->it1, it->it2, cfam);
  }

  return it;
fail:
  it->end = 1;
  return NULL;
}

/* construct an RDL_cycle from the compressed bitset */
RDL_cycle* RDL_cycleIteratorGetCycle(RDL_cycleIterator* it)
{
  unsigned edgeCount = 0, j, mapped_edge;
  RDL_cycle* rdl_cycle = NULL;

  if (it == NULL) {
    RDL_outputFunc(RDL_ERROR, "Iterator is NULL!\n");
    return NULL;
  }

  if (RDL_cycleIteratorAtEnd(it)) {
    RDL_outputFunc(RDL_ERROR, "Cannot retrieve cycle of iterator at end!\n");
    return NULL;
  }

  if (it->mode != 'b') {
    RDL_outputFunc(RDL_ERROR, "Cycle conversion only works for edge defined cycles!\n");
    /* will never happen, when using the public interface! */
    assert(0);
    return NULL;
  }

  rdl_cycle = RDL_alloc(1, sizeof(*rdl_cycle));

  if (!rdl_cycle) {
    goto fail;
  }

  for (j = 0; j < it->data->bccGraphs->bcc_graphs[it->bcc_index]->E; ++j) {
    edgeCount += RDL_bitset_test(it->compressed_cycle, j);
  }
  rdl_cycle->edges = RDL_alloc(edgeCount, sizeof(*rdl_cycle->edges));

  if (!rdl_cycle->edges) {
    goto fail;
  }
  rdl_cycle->urf = it->running_urf;
  rdl_cycle->rcf = it->running_rcf;

  edgeCount = 0;
  for(j=0; j < it->data->bccGraphs->bcc_graphs[it->bcc_index]->E; ++j) {
    if(RDL_bitset_test(it->compressed_cycle, j)) {
      mapped_edge = it->data->bccGraphs->edge_from_bcc_mapping[it->bcc_index][j];
      rdl_cycle->edges[edgeCount][0] = it->data->graph->edges[mapped_edge][0];
      rdl_cycle->edges[edgeCount][1] = it->data->graph->edges[mapped_edge][1];
      ++edgeCount;
    }
  }


  rdl_cycle->weight = edgeCount;

  return rdl_cycle;
fail:
  RDL_deleteCycle(rdl_cycle);
  return NULL;
}

/* check if at end, NULL is at end */
int RDL_cycleIteratorAtEnd(RDL_cycleIterator* it)
{
  if (it == NULL) {
    return 1;
  }
  return it->end;
}

void RDL_deleteCycleIterator(RDL_cycleIterator* it)
{
  if (!it) {
    return;
  }
  if (it->it1) {
    RDL_deletePathIterator(it->it1);
  }
  if (it->it2) {
    RDL_deletePathIterator(it->it2);
  }

  if (it->compressed_cycle) {
    RDL_free(it->compressed_cycle);
  }

  RDL_free(it);
}

/* initialize cycle iterator */
RDL_cycleIterator* RDL_initCycleIterator(
    RDL_IteratorType itype,
    unsigned rcf_index,
    unsigned rcf_index_max,
    unsigned urf_index,
    unsigned urf_index_max,
    unsigned bcc_index,
    unsigned bcc_index_max,
    char mode,
    const RDL_data* data)
{
  unsigned i, j, number;
  RDL_cycleIterator *it;
  const RDL_URFinfo *uInfo;

  it = RDL_alloc(1, sizeof(*it));

  if (!it) {
    goto fail;
  }

  it->mode = mode;
  it->it1 = NULL;
  it->it2 = NULL;
  it->end = 0;
  it->data = data;
  it->compressed_cycle = NULL;

  /* stop right here, if there is nothing to iterate */
  if (it->data->bccGraphs->nof_bcc == 0) {
    it->end = 1;
    return it;
  }
  if (it->mode == 'a') {
    number = it->data->bccGraphs->bcc_graphs[bcc_index]->V;
  }
  else {
    number = it->data->bccGraphs->bcc_graphs[bcc_index]->E;
  }
  if (RDL_bitset_init(&it->compressed_cycle, number) == RDL_INVALID_RESULT) {
    goto fail;
  }

  it->bcc_index = bcc_index;
  it->bcc_index_max = bcc_index_max;

  uInfo = data->urfInfoPerBCC[bcc_index];

  if (itype == RDL_URF_IT || itype == RDL_RCF_IT) {
    it->urf_index = urf_index;
    it->urf_index_max = urf_index_max;
  }
  else {
    it->urf_index = 0;
    it->urf_index_max = data->nofURFsPerBCC[it->bcc_index] - 1;
  }

  if (itype == RDL_RCF_IT) {
    it->rcf_index = rcf_index;
    it->rcf_index_max = rcf_index_max;
  }
  else {
    it->rcf_index = 0;
    it->rcf_index_max = uInfo->nofCFsPerURF[it->urf_index] - 1;
  }

  it->type = itype;
  it->running_rcf = 0;
  it->running_urf = 0;

  for(i = 0; i < bcc_index; ++i) {
    for(j = 0; j < data->nofURFsPerBCC[i]; ++j) {
      it->running_rcf += data->urfInfoPerBCC[i]->nofCFsPerURF[j];
    }
    it->running_urf += data->nofURFsPerBCC[i];
  }
  it->running_urf += urf_index;

  for(i = 0; i < urf_index; ++i) {
    it->running_rcf += data->urfInfoPerBCC[bcc_index]->nofCFsPerURF[i];
  }
  it->running_rcf += rcf_index;

  if (!RDL_cycleIteratorNext(it)) {
    goto fail;
  }

  return it;
fail:
  RDL_deleteCycleIterator(it);
  return NULL;
}

