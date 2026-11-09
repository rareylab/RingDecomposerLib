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

#include "RDLgraph.h"

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>

#include "RDLutility.h"

static void RDL_DFSvisit(const RDL_graph *gra, unsigned vertex, char *visited)
{
  unsigned i,j;
  for(i=0; i<gra->degree[vertex]; ++i)
  {
    j = gra->adjList[vertex][i][0];
    if(visited[j] == 0)
    {
      visited[j] = 1;
      RDL_DFSvisit(gra,j,visited);
    }
  }
}

/** returns 1 if the graph is connected, 0 otherwise. uses DFS */
int RDL_checkGraphConnected(const RDL_graph *gra)
{
  unsigned i;
  char *visited;
  int result;
  result = 1;
  if (!gra->V) {
    return 1;
  }
  visited = RDL_alloc(gra->V, sizeof(*visited));
  if (!visited) {
    return -1;
  }
  visited[0] = 1; /*start vertex*/
  for(i=1; i<gra->V; ++i)
  {
    visited[i] = 0; /*unvisited*/
  }
  RDL_DFSvisit(gra, 0, visited);
  for(i=0; i<gra->V; ++i)
  {
    if(visited[i] == 0) /*one was unvisited*/
    {
      result = 0;
      break;
    }
  }

  RDL_free(visited);
  return result;
}

/* if nodes are adjacent */
int RDL_isAdj(const RDL_graph *graph, unsigned i, unsigned j)
{
  unsigned idx;
  for(idx=0; idx<graph->degree[i]; ++idx)
  {
    if(graph->adjList[i][idx][0] == j)
    {
      return 1;
    }
  }
  return 0;
}

void RDL_deleteGraph(RDL_graph *gra)
{
  unsigned i;
  if (!gra) {
    return;
  }
  if(gra->edges)
  {
    for(i=0; i<gra->E; ++i)
    {
      RDL_free(gra->edges[i]);
    }
    RDL_free(gra->edges);
  }
  for(i=0; i<gra->V; ++i)
  {
    if(gra->adjList)
    {
      RDL_free(gra->adjList[i]);
    }
  }
  RDL_free(gra->adjList);
  RDL_free(gra->degree);
  RDL_free(gra);
}

void RDL_printGraph(const RDL_graph *graph)
{
  unsigned i,j;

  printf("|V|=%d, |E|=%d\n",graph->V,graph->E);
  for(i=0; i<graph->V; ++i)
  {
    printf("%d:  ",i);
    for(j=0; j<graph->degree[i]; ++j)
    {
      printf("%d ",graph->adjList[i][j][0]);
    }
    printf("\n");
  }
  if(graph->edges)/*undirected*/
  {
    printf("edges:\n");
    for(i=0; i<graph->E; ++i)
    {
      printf("%d: [%d,%d]\n", i, graph->edges[i][0], graph->edges[i][1]);
    }
  }
}

RDL_graph *RDL_initNewGraph_g(unsigned V, char owns_edges)
{
  RDL_graph *graph = RDL_alloc(1, sizeof(*graph));
  if (!graph) {
    return NULL;
  }
  graph->V = V;
  graph->degree = RDL_alloc(V, sizeof(*graph->degree));
  graph->adjList = RDL_alloc(V, sizeof(*graph->adjList));
  if (owns_edges) {
    graph->edgesAlloced = 64;
    graph->edges = RDL_alloc(graph->edgesAlloced, sizeof(*graph->edges));
  }
  if (!graph->degree || !graph->adjList || (owns_edges && !graph->edges)) {
    RDL_deleteGraph(graph);
    return NULL;
  }
  return graph;
}

int RDL_addEdge(RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned i;
  RDL_edge *adj;
  if (from == to) {
    return 1;
  }
  for (i = 0; i < gra->degree[from]; ++i) {
    if (gra->adjList[from][i][0] == to) {
      return 1;
    }
  }
  if (gra->degree[from] == UINT_MAX || gra->E == UINT_MAX) {
    return 0;
  }
  adj = RDL_resize(gra->adjList[from], (size_t)gra->degree[from] + 1, sizeof(*adj));
  if (!adj) {
    return 0;
  }
  gra->adjList[from] = adj;
  adj[gra->degree[from]++][0] = to;
  ++gra->E;
  return 1;
}

unsigned RDL_addUEdge_g(RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned i, *edge, **edges;
  RDL_edge *adj;
  if (!gra || from >= gra->V || to >= gra->V || from == to) {
    RDL_outputFunc(RDL_ERROR, "Invalid edge endpoints.\n");
    return RDL_INVALID_RESULT;
  }
  for (i = 0; i < gra->degree[from]; ++i) {
    if (gra->adjList[from][i][0] == to) {
      return RDL_DUPLICATE_EDGE;
    }
  }
  if (gra->E >= UINT_MAX - 3 || gra->degree[from] == UINT_MAX || gra->degree[to] == UINT_MAX) {
    return RDL_INVALID_RESULT;
  }
  /* Reserve everything before changing counts or writing an edge. A failed
   * reserve may grow storage, but leaves the logical graph unchanged. */
  adj = RDL_resize(gra->adjList[from], (size_t)gra->degree[from] + 1, sizeof(*adj));
  if (!adj) {
    return RDL_INVALID_RESULT;
  }
  gra->adjList[from] = adj;
  adj = RDL_resize(gra->adjList[to], (size_t)gra->degree[to] + 1, sizeof(*adj));
  if (!adj) {
    return RDL_INVALID_RESULT;
  }
  gra->adjList[to] = adj;
  if (gra->E == gra->edgesAlloced) {
    if (gra->edgesAlloced > UINT_MAX / 2) {
      return RDL_INVALID_RESULT;
    }
    i = gra->edgesAlloced ? gra->edgesAlloced * 2 : 64;
    edges = RDL_resize(gra->edges, i, sizeof(*edges));
    if (!edges) {
      return RDL_INVALID_RESULT;
    }
    gra->edges = edges;
    gra->edgesAlloced = i;
  }
  edge = RDL_alloc(2, sizeof(*edge));
  if (!edge) {
    return RDL_INVALID_RESULT;
  }
  edge[0] = from < to ? from : to;
  edge[1] = from < to ? to : from;
  i = gra->E++;
  gra->edges[i] = edge;
  gra->adjList[from][gra->degree[from]][0] = to;
  gra->adjList[from][gra->degree[from]++][1] = i;
  gra->adjList[to][gra->degree[to]][0] = from;
  gra->adjList[to][gra->degree[to]++][1] = i;
  return i;
}

unsigned RDL_edgeId(const RDL_graph *gra, unsigned from, unsigned to)
{
  unsigned j, edge;

  if(from > to) {
    /*swap order to make from < to*/
    edge = to;
    to = from;
    from = edge;
  }

  edge = RDL_INVALID_RESULT;

  for(j=0; j<gra->degree[from]; ++j) {
    if(gra->adjList[from][j][0] == to) {
      edge = gra->adjList[from][j][1];
      break;
    }
  }

  return edge;
}
