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

#include "RDLtarjan.h"

#include <assert.h>
#include <limits.h>
#include <stdlib.h>

#include "RDLgraph.h"
#include "RDLstack.h"
#include "RDLutility.h"

typedef struct RDL_dfsStackElement {
    unsigned u;
    unsigned parent;
    unsigned j;
    unsigned expected_time;
} RDL_dfsStackElement;

/*
 * DFS visit of tarjan's algorithm
 *
 * modified to work with an explicit stack:
 * the stack emulates the recursion by visiting nodes
 * twice
 */
static int RDL_tarjanVisit(const RDL_graph* graph, unsigned start, unsigned* d,
    unsigned* low, unsigned* time, unsigned* curr_bcc, unsigned* bcc,
    RDL_stack* edge_stack)
{
  unsigned j, v, curr_edge, uv_edge, u, parent;
  unsigned *element;
  unsigned *edge_elements = NULL;
  unsigned next_free_edge_element = 0;

  RDL_dfsStackElement *dfs_element, *dfs_element_new;
  RDL_stack* dfs_stack;

  RDL_dfsStackElement *dfs_elements = NULL;
  unsigned next_free_dfs_element = 0;

  dfs_stack = RDL_stack_new();
  if (!dfs_stack) {
    goto fail;
  }

  dfs_elements = RDL_alloc(graph->V, sizeof(*dfs_elements));
  if (!dfs_elements) {
    goto fail;
  }
  dfs_element = &(dfs_elements[next_free_dfs_element]);
  ++next_free_dfs_element;

  dfs_element->u = start;
  dfs_element->parent = UINT_MAX;
  dfs_element->j = 0;
  dfs_element->expected_time = UINT_MAX;
  ++(*time);
  d[start] = *time;
  low[start] = d[start];

  if (!RDL_stack_push(dfs_stack, dfs_element)) {
    goto fail;
  }

  edge_elements = RDL_alloc(graph->E, sizeof(*edge_elements));
  if (!edge_elements) {
    goto fail;
  }

  while (!RDL_stack_empty(dfs_stack)) {

    dfs_element = RDL_stack_top(dfs_stack);
    u = dfs_element->u;
    parent = dfs_element->parent;
    j = dfs_element->j;

    /* while there is one node left in it's adjacency list */
    if (j < graph->degree[u]) {
      v = graph->adjList[u][j][0];
      uv_edge = RDL_edgeId(graph, u, v);

      /* if looking at this node this neighbour the first time */
      if (d[v] == 0) {
        /* DO NOT ADVANCE j, we'll be looking at this node again */
        /* push edge onto edge_stack as in tarjan's */
        edge_elements[next_free_edge_element] = uv_edge;
        if (!RDL_stack_push(edge_stack, &(edge_elements[next_free_edge_element]))) {
          goto fail;
        }
        ++next_free_edge_element;

        ++(*time);
        d[v] = *time;
        low[v] = d[v];

        /* and recurse */
        dfs_element_new = &(dfs_elements[next_free_dfs_element]);
        ++next_free_dfs_element;
        dfs_element_new->u = v;
        dfs_element_new->parent = u;
        dfs_element_new->j = 0;
        dfs_element_new->expected_time = UINT_MAX;
        if (!RDL_stack_push(dfs_stack, dfs_element_new)) {
          goto fail;
        }
        /*
         * THIS IS IMPORTANT: this is the time we expect for v to
         * be discovered when we come back to u
         */
        dfs_element->expected_time = (*time);
      }
      else {
        /* in all other cases, we'll NOT be looking again at this neighbour */
        ++dfs_element->j;

        /*
         * this is interesting: if we're looking at this neighbour the
         * second time and d[v] == d[u] + 1 we're back from recursion
         * and this neighbour was discovered when recursing from u as parent
         */
        if (d[v] == dfs_element->expected_time) {
          /* do tarjan's stuff */
          low[u] = low[u] < low[v] ? low[u] : low[v];
          if (low[v] >= d[u]) {
            do {
              element = RDL_stack_top(edge_stack);
              curr_edge = *element;
              RDL_stack_pop(edge_stack);

              bcc[curr_edge] = *curr_bcc;
            } while (curr_edge != uv_edge);
            ++(*curr_bcc);
          }
        }
        /*
         * this is looking at v as u's neighbour the first time, but
         * v was already discovered
         */
        else if (d[v] < d[u] && v != parent) {
          /* do tarjan's stuff */
          edge_elements[next_free_edge_element] = uv_edge;
          if (!RDL_stack_push(edge_stack, &(edge_elements[next_free_edge_element]))) {
            goto fail;
          }
          ++next_free_edge_element;
          low[u] = low[u] < d[v] ? low[u] : d[v];
        }
      }
    }
    else {
      /* if there is no neighbour left => u is finished and goes from stack */
      RDL_stack_pop(dfs_stack);
    }
  }

  RDL_free(edge_elements);
  RDL_free(dfs_elements);
  RDL_stack_delete(dfs_stack);
  return 1;
fail:
  RDL_free(edge_elements);
  RDL_free(dfs_elements);
  RDL_stack_delete(dfs_stack);
  return 0;
}

/* perform tarjan's algorithm */
static unsigned RDL_tarjan(const RDL_graph* graph,
                    unsigned** bcc)
{
  unsigned *d = NULL, *low = NULL, u, time, curr_bcc, i;
  unsigned result = RDL_INVALID_RESULT;
  RDL_stack* edge_stack = NULL;

  *bcc = NULL;
  d = RDL_alloc(graph->V, sizeof(*d));
  if (!d) {
    goto fail;
  }
  low = RDL_alloc(graph->V, sizeof(*low));
  if (!low) {
    goto fail;
  }
  *bcc = RDL_alloc(graph->E, sizeof(**bcc));
  if (!*bcc) {
    goto fail;
  }
  edge_stack = RDL_stack_new();
  if (!edge_stack) {
    goto fail;
  }

  time = 0;
  curr_bcc = 1;

  for (u = 0; u < graph->V; ++u) {
    if (d[u] == 0) {
      if (!RDL_tarjanVisit(graph, u, d, low, &time, &curr_bcc, *bcc, edge_stack)) {
        goto fail;
      }
    }
  }

  result = curr_bcc - 1;
fail:
  RDL_free(d);
  RDL_free(low);
  RDL_stack_delete(edge_stack);
  if (result == RDL_INVALID_RESULT) {
    RDL_free(*bcc);
    *bcc = NULL;
  }
  return result;
}

/* perform Tarjan's algorithm and save mapping */
RDL_BCCGraph* RDL_tarjanBCC(const RDL_graph* graph)
{
  unsigned *bcc = NULL, u, node, i, j, k, ii[2],
    nof_bcc, *bcc_size = NULL, nof_nontrivial_bcc=0, curr_bcc, found;
  unsigned *non_trivial_mapping = NULL;
  void *resized;
  RDL_BCCGraph* result = NULL;

  nof_bcc = RDL_tarjan(graph, &bcc);
  if (nof_bcc == RDL_INVALID_RESULT) {
    return NULL;
  }

  result = RDL_alloc(1, sizeof(*result));
  if (!result) {
    goto fail;
  }
  result->complete_graph = graph;
  bcc_size = RDL_alloc(nof_bcc, sizeof(*bcc_size));
  if (!bcc_size) {
    goto fail;
  }
  non_trivial_mapping = RDL_alloc(nof_bcc, sizeof(*non_trivial_mapping));
  if (!non_trivial_mapping) {
    goto fail;
  }

  for (i = 0; i < nof_bcc; ++i) {
    bcc_size[i] = 0;
  }

  for (i = 0; i < graph->E; ++i) {
    if (bcc[i] == 0) {
      RDL_outputFunc(RDL_ERROR, "edge has no associated BCC %u\n", i);
      assert(0); /* should NEVER happen */
    }
    curr_bcc = bcc[i];
    ++bcc_size[curr_bcc - 1];
  }

  for (i = 0; i < nof_bcc; ++i) {
    if(bcc_size[i] > 1) {
      non_trivial_mapping[i] = nof_nontrivial_bcc;
      ++nof_nontrivial_bcc;
    }
    else {
      non_trivial_mapping[i] = RDL_NO_RINGSYSTEM;
    }
  }

  result->nof_bcc = nof_nontrivial_bcc;
  result->bcc_graphs = RDL_alloc(nof_nontrivial_bcc, sizeof(*result->bcc_graphs));
  if (!result->bcc_graphs) {
    goto fail;
  }
  result->edge_to_bcc_mapping = RDL_alloc(graph->E, sizeof(*result->edge_to_bcc_mapping));
  if (!result->edge_to_bcc_mapping) {
    goto fail;
  }

  for (i = 0; i < graph->E; ++i) {
    result->edge_to_bcc_mapping[i] = RDL_alloc(2, sizeof(**result->edge_to_bcc_mapping));
    if (!result->edge_to_bcc_mapping[i]) {
      goto fail;
    }
    result->edge_to_bcc_mapping[i][0] = RDL_NO_RINGSYSTEM;
    result->edge_to_bcc_mapping[i][1] = RDL_NO_RINGSYSTEM;
  }

  result->edge_from_bcc_mapping = RDL_alloc(nof_nontrivial_bcc, sizeof(*result->edge_from_bcc_mapping));
  if (!result->edge_from_bcc_mapping) {
    goto fail;
  }
  result->nof_edges_per_bcc = RDL_alloc(nof_nontrivial_bcc, sizeof(*result->nof_edges_per_bcc));
  if (!result->nof_edges_per_bcc) {
    goto fail;
  }
  for (i = 0; i < nof_nontrivial_bcc; ++i) {
    result->edge_from_bcc_mapping[i] = NULL;
    result->nof_edges_per_bcc[i] = 0;
  }

  result->node_to_bcc_mapping = RDL_alloc(graph->V, sizeof(*result->node_to_bcc_mapping));
  if (!result->node_to_bcc_mapping) {
    goto fail;
  }
  result->nof_bcc_per_node = RDL_alloc(graph->V, sizeof(*result->nof_bcc_per_node));
  if (!result->nof_bcc_per_node) {
    goto fail;
  }
  for (i = 0; i < graph->V; ++i) {
    result->nof_bcc_per_node[i] = 0;
    result->node_to_bcc_mapping[i] = NULL;
  }

  result->node_from_bcc_mapping = RDL_alloc(nof_nontrivial_bcc, sizeof(*result->node_from_bcc_mapping));
  if (!result->node_from_bcc_mapping) {
    goto fail;
  }
  result->nof_nodes_per_bcc = RDL_alloc(nof_nontrivial_bcc, sizeof(*result->nof_nodes_per_bcc));
  if (!result->nof_nodes_per_bcc) {
    goto fail;
  }
  for (i = 0; i < nof_nontrivial_bcc; ++i) {
    result->node_from_bcc_mapping[i] = NULL;
    result->nof_nodes_per_bcc[i] = 0;
  }

  for (i = 0; i < graph->E; ++i) {
    curr_bcc = bcc[i];

    /* skip trivial BCCs */
    if (bcc_size[curr_bcc-1] <= 1) {
      continue;
    }

    /* -1 because of 1-based BCCs */
    u = non_trivial_mapping[curr_bcc - 1];
    result->edge_to_bcc_mapping[i][0] = u;
    result->edge_to_bcc_mapping[i][1] = result->nof_edges_per_bcc[u]++;
    resized = RDL_resize(result->edge_from_bcc_mapping[u], result->nof_edges_per_bcc[u], sizeof(**result->edge_from_bcc_mapping));
    if (!resized) {
      goto fail;
    }
    result->edge_from_bcc_mapping[u] = resized;
    result->edge_from_bcc_mapping[u][result->nof_edges_per_bcc[u] - 1] = i;

    for (k = 0; k <= 1; ++k) {
      node = graph->edges[i][k];

      found = 0;
      for (j = 0; j < result->nof_bcc_per_node[node] && !found; ++j) {
        if (result->node_to_bcc_mapping[node][2 * j] == u) {
          found = 1;
        }
      }

      if (!found) {
        ++result->nof_bcc_per_node[node];
        resized = RDL_resize(result->node_to_bcc_mapping[node], (size_t)2 * result->nof_bcc_per_node[node], sizeof(*result->node_to_bcc_mapping[node]));
        if (!resized) {
          goto fail;
        }
        result->node_to_bcc_mapping[node] = resized;
        result->node_to_bcc_mapping[node][2*(result->nof_bcc_per_node[node] - 1)] = u;
        result->node_to_bcc_mapping[node][2*(result->nof_bcc_per_node[node] - 1)+1] = result->nof_nodes_per_bcc[u]++;

        resized = RDL_resize(result->node_from_bcc_mapping[u], result->nof_nodes_per_bcc[u], sizeof(**result->node_from_bcc_mapping));
        if (!resized) {
          goto fail;
        }
        result->node_from_bcc_mapping[u] = resized;
        result->node_from_bcc_mapping[u][result->nof_nodes_per_bcc[u]-1] = node;
      }
    }
  }

  for (u = 0; u < nof_nontrivial_bcc; ++u) {
    result->bcc_graphs[u] = RDL_initNewGraph(result->nof_nodes_per_bcc[u]);
    if (!result->bcc_graphs[u]) {
      goto fail;
    }
  }

  for (i = 0; i < graph->E; ++i) {
    u = result->edge_to_bcc_mapping[i][0];

    /* skip if there is no BCC */
    if (u == RDL_NO_RINGSYSTEM) {
      continue;
    }

    ii[0] = ii[1] = RDL_NO_RINGSYSTEM;

    for (k = 0; k <= 1; ++ k) {
      for (j = 0; j < result->nof_bcc_per_node[graph->edges[i][k]]; ++j) {
        if (result->node_to_bcc_mapping[graph->edges[i][k]][2 * j] == u) {
          ii[k] = result->node_to_bcc_mapping[graph->edges[i][k]][2 * j + 1];
        }
      }
      if (ii[k] == RDL_NO_RINGSYSTEM) {
        RDL_outputFunc(RDL_ERROR, "node %u not part of BCC %u!\n",
            graph->edges[i][k], u);
        /* should never happen */
        assert(0);
      }
    }
    if (RDL_addUEdge(result->bcc_graphs[u], ii[0], ii[1]) == RDL_INVALID_RESULT) {
      goto fail;
    }
  }

  RDL_free(bcc);
  RDL_free(bcc_size);
  RDL_free(non_trivial_mapping);

  return result;
fail:
  RDL_free(bcc); RDL_free(bcc_size); RDL_free(non_trivial_mapping);
  RDL_deleteBCCGraph(result);
  return NULL;
}

void RDL_deleteBCCGraph(RDL_BCCGraph* graph)
{
  unsigned i;

  if (!graph) {
    return;
  }
  for (i = 0; i < graph->nof_bcc; ++i) {
    if (graph->bcc_graphs) {
      RDL_deleteGraph(graph->bcc_graphs[i]);
    }
    if (graph->edge_from_bcc_mapping) {
      RDL_free(graph->edge_from_bcc_mapping[i]);
    }
    if (graph->node_from_bcc_mapping) {
      RDL_free(graph->node_from_bcc_mapping[i]);
    }
  }
  RDL_free(graph->bcc_graphs);
  RDL_free(graph->edge_from_bcc_mapping);
  RDL_free(graph->node_from_bcc_mapping);
  RDL_free(graph->nof_edges_per_bcc);
  RDL_free(graph->nof_nodes_per_bcc);

  for (i = 0; i < graph->complete_graph->V; ++i) {
    if (graph->node_to_bcc_mapping) {
      RDL_free(graph->node_to_bcc_mapping[i]);
    }
  }
  RDL_free(graph->node_to_bcc_mapping);

  for (i = 0; i < graph->complete_graph->E; ++i) {
    if (graph->edge_to_bcc_mapping) {
      RDL_free(graph->edge_to_bcc_mapping[i]);
    }
  }
  RDL_free(graph->edge_to_bcc_mapping);
  RDL_free(graph->nof_bcc_per_node);
  RDL_free(graph);
}
