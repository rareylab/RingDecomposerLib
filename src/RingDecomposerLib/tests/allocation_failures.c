/* Fail each allocation in turn and verify ownership, error returns and cleanup. */
#include <stdlib.h>
#include <stdio.h>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <assert.h>
#include <string.h>
#include "failure_alloc.h"
#include "RingDecomposerLib.h"
#include "RDLgraph.h"
#include "RDLdimacs.h"
#include "RDLtesting.h"
#include "RDLutility.h"
#include "RDLstack.h"
#undef malloc
#undef calloc
#undef realloc
#undef free

static long fail_at = -1, calls, live;
static unsigned errors;
static void report_error(RDL_ERROR_LEVEL level, const char *fmt, ...)
{
  (void)fmt;
  if (level == RDL_ERROR) {
    ++errors;
  }
}
static void *allocations[100000];
static void record(void *p)
{
  unsigned i;
  if (!p) {
    return;
  }
  for (i = 0; i < 100000; ++i) {
    if (!allocations[i]) {
      allocations[i] = p;
      ++live;
      return;
    }
  }
  abort();
}
static void forget(void *p)
{
  unsigned i;
  if (!p) {
    return;
  }
  for (i = 0; i < 100000; ++i) {
    if (allocations[i] == p) {
      allocations[i] = NULL;
      --live;
      return;
    }
  }
  fprintf(stderr, "Attempt to free unowned memory\n");
  abort();
}
void *audit_malloc(size_t n)
{
  void *p;
  if (calls++ == fail_at) {
    return NULL;
  }
  p = malloc(n); record(p); return p;
}
void *audit_calloc(size_t n, size_t size)
{
  void *p;
  if (calls++ == fail_at) {
    return NULL;
  }
  p = calloc(n, size); record(p); return p;
}
void *audit_realloc(void *p, size_t n)
{
  void *q;
  if (calls++ == fail_at) {
    return NULL;
  }
  /* Forget before realloc: p may no longer be valid afterwards. */
  forget(p);
  q = realloc(p, n);
  record(q ? q : p);
  return q;
}
void audit_free(void *p) { forget(p); free(p); }

static RDL_graph *make_graph(unsigned shape)
{
  unsigned i, j;
  RDL_graph *g = RDL_initNewGraph(shape == 3 ? 75 : (shape == 5 ? 9 : 8));
  assert(g);
  if (shape == 0) {
    return g;
  }
  if (shape == 1 || shape == 5) {
    for (i = 0; i < (shape == 5 ? 9 : 8); ++i) {
      for (j = i + 1; j < (shape == 5 ? 9 : 8); ++j) {
        assert(RDL_addUEdge(g, i, j) != RDL_INVALID_RESULT);
      }
    }
  } else if (shape == 2) {
    for (i = 0; i < 6; i += 3) {
      RDL_addUEdge(g, i, i+1); RDL_addUEdge(g, i+1, i+2);
      RDL_addUEdge(g, i+2, i);
    }
  } else if (shape == 3) {
    for (i = 0; i < 75; ++i) {
      RDL_addUEdge(g, i, (i + 1) % 75);
    }
  }
  if (shape == 4) {
    for (i = 0; i < 3; ++i) {
      for (j = 3; j < 6; ++j) {
        RDL_addUEdge(g, i, j);
      }
    }
    RDL_addUEdge(g, 5, 6); RDL_addUEdge(g, 6, 7); RDL_addUEdge(g, 7, 3);
  }
  return g;
}

static int query(RDL_data *data, unsigned op)
{
  unsigned n = 0;
  RDL_cycle **cycles = NULL;
  RDL_edge *edges = NULL;
  RDL_node *nodes = NULL;
  unsigned *families = NULL;
  char **bits = NULL;
  double count;
  RDL_cycleIterator *it;
  switch (op) {
    case 0: n = RDL_getSSSR(data, &cycles); break;
    case 1: n = RDL_getRCPrototypes(data, &cycles); break;
    case 2: n = RDL_getRCycles(data, &cycles); break;
    case 3: n = RDL_getNodesForURF(data, 0, &nodes); break;
    case 4: n = RDL_getNodesForRCF(data, 0, &nodes); break;
    case 5: n = RDL_getEdgesForURF(data, 0, &edges); break;
    case 6: n = RDL_getEdgesForRCF(data, 0, &edges); break;
    case 7: n = RDL_getURFsContainingNode(data, 0, &families); break;
    case 8: n = RDL_getRCFsContainingNode(data, 0, &families); break;
    case 9: n = RDL_getURFsContainingEdge(data, 0, 1, &families); break;
    case 10: n = RDL_getRCFsContainingEdge(data, 0, 1, &families); break;
    case 11: n = RDL_getEdgeArray(data, &edges); break;
    case 12: n = RDL_getEdgesForRingsystem(data, 0, &edges); break;
    case 13: n = RDL_getNodesForRingsystem(data, 0, &nodes); break;
    case 14:
      count = RDL_getNofRC(data);
      return count != RDL_INVALID_RC_COUNT;
    case 15:
      it = RDL_getRCyclesIterator(data);
      if (!it) {
        return 0;
      }
      while (!RDL_cycleIteratorAtEnd(it)) {
        if (!RDL_cycleIteratorNext(it)) { RDL_deleteCycleIterator(it); return 0; }
      }
      RDL_deleteCycleIterator(it);
      return 1;
    case 16:
      n = RDL_getRCPrototypes(data, &cycles);
      if (n != RDL_INVALID_RESULT) {
        unsigned translated = RDL_translateCycArray(data, cycles, n, &bits);
        RDL_deleteCycles(cycles, n);
        if (translated == RDL_INVALID_RESULT) {
          return 0;
        }
        RDL_delete2DCharArray(bits, translated);
        return 1;
      }
      break;
    case 17: return RDL_checkConsistency(data) == 0;
    default: abort();
  }
  if (n == RDL_INVALID_RESULT) {
    assert(!cycles && !edges && !nodes && !families);
    return 0;
  }
  RDL_deleteCycles(cycles, cycles ? n : 0);
  audit_free(edges); audit_free(nodes); audit_free(families);
  return 1;
}

static void sweep_calculation(unsigned shape)
{
  long point, total, baseline;
  for (point = 0; ; ++point) {
    RDL_graph *g = make_graph(shape);
    RDL_data *data;
    baseline = live; calls = 0; errors = 0; fail_at = point;
    data = RDL_calculate(g);
    total = calls; fail_at = -1;
    if (point < total) { assert(!data); assert(errors); }
    if (data) {
      RDL_deleteData(data);
    } else {
      assert(live == baseline);
      RDL_deleteGraph(g);
    }
    assert(live == 0);
    if (point >= total) {
      break;
    }
  }
  printf("Calculation shape %u: %ld allocation points\n", shape, total);
}

static void sweep_queries(unsigned shape)
{
  unsigned op;
  for (op = 0; op < (shape == 2 ? 18 : 17); ++op) {
    long point, total, baseline;
    for (point = 0; ; ++point) {
      RDL_data *data = RDL_calculate(make_graph(shape));
      int success;
      assert(data);
      baseline = live; calls = 0; errors = 0; fail_at = point;
      success = query(data, op);
      total = calls; fail_at = -1;
      if (point < total) {
        assert(errors);
      }
      if (point < total && success) {
        fprintf(stderr, "Query %u shape %u swallowed failure %ld\n", op, shape, point);
        abort();
      }
      if (live != baseline) {
        fprintf(stderr, "Query %u shape %u failure %ld: live %ld, expected %ld\n",
                op, shape, point, live, baseline); abort();
      }
      RDL_deleteData(data);
      assert(live == 0);
      if (point >= total) {
        break;
      }
    }
    printf("Query %u shape %u: %ld allocation points\n", op, shape, total);
  }
}

static void sweep_validators(void)
{
  const char a[] = {1, 1, 0, 0}, b[] = {0, 1, 1, 0}, c[] = {1, 0, 1, 0};
  const char *rings[] = {a, b, c};
  unsigned bcc[] = {0, 0, 0}, families[] = {0, 0, 0}, op;
  for (op = 0; op < 2; ++op) {
    long point, total;
    for (point = 0; ; ++point) {
      int result;
      calls = 0; errors = 0; fail_at = point;
      if (op == 0) {
        result = RDL_validateBasis(rings, 3, rings, 2, 4, bcc, RDL_no_stop_fun);
      } else {
        result = RDL_validateRingFamilies(rings, families, 3, 4, RDL_no_stop_fun);
      }
      total = calls; fail_at = -1;
      if (point < total) {
        assert(result == -2);
      }
      assert(live == 0);
      if (point >= total) {
        break;
      }
    }
    printf("Validator %u: %ld allocation points\n", op, total);
  }
}

static void sweep_insertion(void)
{
  unsigned full;
  for (full = 0; full < 2; ++full) {
    long point, total;
    for (point = 0; ; ++point) {
      RDL_graph *g = RDL_initNewGraph(70);
      unsigned i, result, before;
      if (full) {
        for (i = 0; i < 64; ++i) {
          RDL_addUEdge(g, i, i + 1);
        }
      }
      before = g->E;
      calls = 0; errors = 0; fail_at = point;
      result = RDL_addUEdge(g, 65, 66);
      total = calls; fail_at = -1;
      if (point < total) {
        assert(result == RDL_INVALID_RESULT && g->E == before);
        assert(g->degree[65] == 0 && g->degree[66] == 0);
        assert(RDL_addUEdge(g, 65, 66) != RDL_INVALID_RESULT);
      }
      RDL_deleteGraph(g); assert(live == 0);
      if (point >= total) {
        break;
      }
    }
  }
}

int main(int argc, char **argv)
{
  unsigned shape;
  long point, total;
  assert(argc == 2);
  RDL_setOutputFunction(report_error);
  {
    void *original = RDL_alloc(2, sizeof(unsigned));
    assert(original);
    assert(!RDL_resize(original, (size_t)-1, 2));
    assert(errors);
    audit_free(original);
    assert(!RDL_alloc((size_t)-1, 2));
  }
  for (point = 0; ; ++point) {
    RDL_graph *g;
    calls = 0; errors = 0; fail_at = point;
    g = RDL_initNewGraph(3);
    total = calls; fail_at = -1;
    if (point < total) { assert(!g); assert(errors); }
    RDL_deleteGraph(g); assert(live == 0);
    if (point >= total) {
      break;
    }
  }
  sweep_validators();
  sweep_insertion();
  for (shape = 0; shape < 6; ++shape) {
    sweep_calculation(shape);
  }
  for (shape = 1; shape < 5; ++shape) {
    sweep_queries(shape);
  }
  for (point = 0; ; ++point) {
    RDL_DimacsGraph *g;
    calls = 0; errors = 0; fail_at = point;
    g = RDL_dimacsGraph_read(argv[1]);
    total = calls; fail_at = -1;
    if (point < total) { assert(!g); assert(errors); }
    RDL_dimacsGraph_delete(g); assert(live == 0);
    if (point >= total) {
      break;
    }
  }
  printf("DIMACS: %ld allocation points\n", total);
  return 0;
}
