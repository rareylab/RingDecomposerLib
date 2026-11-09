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

enum GraphShape {
  GRAPH_EMPTY,
  GRAPH_COMPLETE_8,
  GRAPH_TWO_TRIANGLES,
  GRAPH_CYCLE_75,
  GRAPH_MIXED,
  GRAPH_COMPLETE_9,
  GRAPH_SHAPE_COUNT
};

static void add_edge_checked(RDL_graph *graph, unsigned from, unsigned to)
{
  unsigned result = RDL_addUEdge(graph, from, to);
  assert(result != RDL_INVALID_RESULT);
}

static RDL_graph *make_graph(unsigned shape)
{
  unsigned i, j, vertex_count;
  RDL_graph *g;

  switch (shape) {
    case GRAPH_EMPTY:
    case GRAPH_COMPLETE_8:
    case GRAPH_TWO_TRIANGLES:
    case GRAPH_MIXED:
      vertex_count = 8;
      break;
    case GRAPH_CYCLE_75:
      vertex_count = 75;
      break;
    case GRAPH_COMPLETE_9:
      vertex_count = 9;
      break;
    default:
      abort();
  }
  g = RDL_initNewGraph(vertex_count);
  assert(g);
  if (shape == GRAPH_EMPTY) {
    return g;
  }
  if (shape == GRAPH_COMPLETE_8 || shape == GRAPH_COMPLETE_9) {
    for (i = 0; i < vertex_count; ++i) {
      for (j = i + 1; j < vertex_count; ++j) {
        add_edge_checked(g, i, j);
      }
    }
  } else if (shape == GRAPH_TWO_TRIANGLES) {
    for (i = 0; i < 6; i += 3) {
      add_edge_checked(g, i, i+1);
      add_edge_checked(g, i+1, i+2);
      add_edge_checked(g, i+2, i);
    }
  } else if (shape == GRAPH_CYCLE_75) {
    for (i = 0; i < vertex_count; ++i) {
      add_edge_checked(g, i, (i + 1) % vertex_count);
    }
  }
  if (shape == GRAPH_MIXED) {
    for (i = 0; i < 3; ++i) {
      for (j = 3; j < 6; ++j) {
        add_edge_checked(g, i, j);
      }
    }
    add_edge_checked(g, 5, 6);
    add_edge_checked(g, 6, 7);
    add_edge_checked(g, 7, 3);
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

typedef int (*SweepOperation)(void *context, long point);

static void sweep(const char *name, SweepOperation operation, void *context)
{
  long point, total;
  for (point = 0; ; ++point) {
    calls = 0;
    errors = 0;
    fail_at = -1;
    (void)operation(context, point);
    total = calls;
    fail_at = -1;
    if (point >= total) {
      break;
    }
  }
  printf("%s: %ld allocation points\n", name, total);
}

static int calculation_step(void *context, long point)
{
  unsigned shape = *(unsigned *)context;
  long baseline;
  RDL_graph *g;
  RDL_data *data;
  fail_at = -1;
  g = make_graph(shape);
  baseline = live;
  calls = 0; errors = 0; fail_at = point;
  data = RDL_calculate(g);
  if (point < calls) { assert(!data); assert(errors); }
  if (data) {
    RDL_deleteData(data);
  } else {
    assert(live == baseline);
    RDL_deleteGraph(g);
  }
  fail_at = -1;
  assert(live == 0);
  return data != NULL;
}

typedef struct { unsigned shape, op; } QueryContext;
static int query_step(void *context, long point)
{
  QueryContext *ctx = context;
  long baseline, total;
  RDL_data *data;
  int success;
  fail_at = -1;
  data = RDL_calculate(make_graph(ctx->shape));
  assert(data);
  baseline = live; calls = 0; errors = 0; fail_at = point;
  success = query(data, ctx->op);
  total = calls;
  fail_at = -1;
  if (point < total) {
    assert(errors);
    if (success) {
      fprintf(stderr, "Query %u shape %u swallowed failure %ld\n",
              ctx->op, ctx->shape, point);
      abort();
    }
  }
  if (live != baseline) {
    fprintf(stderr, "Query %u shape %u failure %ld: live %ld, expected %ld\n",
            ctx->op, ctx->shape, point, live, baseline);
    abort();
  }
  RDL_deleteData(data);
  assert(live == 0);
  return success;
}

typedef struct {
  const char **rings;
  unsigned *bcc;
  unsigned *families;
  unsigned op;
} ValidatorContext;
static int validator_step(void *context, long point)
{
  ValidatorContext *ctx = context;
  int result;
  calls = 0; errors = 0; fail_at = point;
  if (ctx->op == 0) {
    result = RDL_validateBasis(ctx->rings, 3, ctx->rings, 2, 4,
        ctx->bcc, RDL_no_stop_fun);
  } else {
    result = RDL_validateRingFamilies(ctx->rings, ctx->families, 3, 4,
        RDL_no_stop_fun);
  }
  if (point < calls) {
    assert(result == -2);
  }
  assert(live == 0);
  return result != -2;
}

static void sweep_validators(void)
{
  const char a[] = {1, 1, 0, 0}, b[] = {0, 1, 1, 0}, c[] = {1, 0, 1, 0};
  const char *rings[] = {a, b, c};
  unsigned bcc[] = {0, 0, 0}, families[] = {0, 0, 0}, op;
  ValidatorContext context;
  char name[40];
  context.rings = rings; context.bcc = bcc; context.families = families;
  for (op = 0; op < 2; ++op) {
    context.op = op;
    snprintf(name, sizeof(name), "Validator %u", op);
    sweep(name, validator_step, &context);
  }
}

typedef struct { unsigned full; } InsertionContext;
static int insertion_step(void *context, long point)
{
  InsertionContext *ctx = context;
  RDL_graph *g;
  unsigned i, result, before;
  fail_at = -1;
  g = RDL_initNewGraph(70);
  assert(g);
  if (ctx->full) {
    for (i = 0; i < 64; ++i) {
      assert(RDL_addUEdge(g, i, i + 1) != RDL_INVALID_RESULT);
    }
  }
  before = g->E;
  calls = 0; errors = 0; fail_at = point;
  result = RDL_addUEdge(g, 65, 66);
  if (point < calls) {
    assert(result == RDL_INVALID_RESULT && g->E == before);
    assert(g->degree[65] == 0 && g->degree[66] == 0);
    fail_at = -1;
    assert(RDL_addUEdge(g, 65, 66) != RDL_INVALID_RESULT);
  }
  fail_at = -1;
  RDL_deleteGraph(g);
  assert(live == 0);
  return result != RDL_INVALID_RESULT;
}

static int graph_init_step(void *context, long point)
{
  RDL_graph *g;
  (void)context;
  calls = 0; errors = 0; fail_at = point;
  g = RDL_initNewGraph(3);
  if (point < calls) { assert(!g); assert(errors); }
  RDL_deleteGraph(g);
  fail_at = -1;
  assert(live == 0);
  return g != NULL;
}

static int dimacs_step(void *context, long point)
{
  const char *path = context;
  RDL_DimacsGraph *g;
  calls = 0; errors = 0; fail_at = point;
  g = RDL_dimacsGraph_read(path);
  if (point < calls) { assert(!g); assert(errors); }
  RDL_dimacsGraph_delete(g);
  fail_at = -1;
  assert(live == 0);
  return g != NULL;
}

int main(int argc, char **argv)
{
  unsigned shape;
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
  sweep("Graph initialization", graph_init_step, NULL);
  sweep_validators();
  {
    InsertionContext context;
    context.full = 0;
    sweep("Edge insertion", insertion_step, &context);
    context.full = 1;
    sweep("Edge insertion with capacity growth", insertion_step, &context);
  }
  for (shape = GRAPH_EMPTY; shape < GRAPH_SHAPE_COUNT; ++shape) {
    char name[64];
    snprintf(name, sizeof(name), "Calculation shape %u", shape);
    sweep(name, calculation_step, &shape);
  }
  for (shape = GRAPH_COMPLETE_8; shape <= GRAPH_MIXED; ++shape) {
    QueryContext context;
    unsigned op;
    context.shape = shape;
    for (op = 0; op < (shape == GRAPH_TWO_TRIANGLES ? 18 : 17); ++op) {
      char name[80];
      context.op = op;
      snprintf(name, sizeof(name), "Query %u shape %u", op, shape);
      sweep(name, query_step, &context);
    }
  }
  sweep("DIMACS", dimacs_step, argv[1]);
  return 0;
}
