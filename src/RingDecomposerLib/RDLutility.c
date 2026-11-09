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

#include "RDLutility.h"

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef RDL_TEST_ALLOCATORS
#include "tests/failure_alloc.h"
#endif

RDL_outputFunction RDL_outputFunc = RDL_writeNothing;

/* Zero-length arrays still have a freeable, non-NULL representation. */
void *RDL_alloc(size_t count, size_t size)
{
  void *result;
  if (size && count > (size_t)-1 / size) {
    RDL_outputFunc(RDL_ERROR, "Allocation size overflow.\n");
    return NULL;
  }
#ifdef RDL_TEST_ALLOCATORS
  result = audit_calloc(count ? count : 1, size ? size : 1);
#else
  result = calloc(count ? count : 1, size ? size : 1);
#endif
  if (!result) {
    RDL_outputFunc(RDL_ERROR, "Memory allocation failed.\n");
  }
  return result;
}

void *RDL_resize(void *ptr, size_t count, size_t size)
{
  void *result;
  if (size && count > (size_t)-1 / size) {
    RDL_outputFunc(RDL_ERROR, "Allocation size overflow.\n");
    return NULL;
  }
#ifdef RDL_TEST_ALLOCATORS
  result = audit_realloc(ptr, (count ? count : 1) * (size ? size : 1));
#else
  result = realloc(ptr, (count ? count : 1) * (size ? size : 1));
#endif
  if (!result) {
    RDL_outputFunc(RDL_ERROR, "Memory allocation failed.\n");
  }
  return result;
}

void RDL_free(void *ptr)
{
#ifdef RDL_TEST_ALLOCATORS
  audit_free(ptr);
#else
  free(ptr);
#endif
}

unsigned **RDL_alloc2DUIntArray(unsigned n, unsigned m)
{
  unsigned i;
  unsigned **arr = RDL_alloc(n, sizeof(*arr));
  if (!arr) {
    return NULL;
  }
  for(i=0; i<n; ++i)
  {
    arr[i] = RDL_alloc(m, sizeof(**arr));
    if (!arr[i]) {
      while (i) {
        RDL_free(arr[--i]);
      }
      RDL_free(arr);
      return NULL;
    }
  }
  return arr;
}

char **RDL_alloc2DCharArray(unsigned n, unsigned m)
{
  unsigned i;
  char **arr = RDL_alloc(n, sizeof(*arr));
  if (!arr) {
    return NULL;
  }
  for(i=0; i<n; ++i)
  {
    arr[i] = RDL_alloc(m, sizeof(**arr));
    if (!arr[i]) {
      while (i) {
        RDL_free(arr[--i]);
      }
      RDL_free(arr);
      return NULL;
    }
  }
  return arr;
}

void RDL_delete2DUIntArray(unsigned **arr, unsigned n)
{
  unsigned i;
  if (!arr) {
    return;
  }
  for (i = 0; i < n; ++i) {
    RDL_free(arr[i]);
  }
  RDL_free(arr);
}

void RDL_delete2DCharArray(char **arr, unsigned n)
{
  unsigned i;
  if (!arr) {
    return;
  }
  for (i = 0; i < n; ++i) {
    RDL_free(arr[i]);
  }
  RDL_free(arr);
}

void RDL_print2DUIntArray(unsigned **arr, unsigned n, unsigned m)
{
  unsigned i,j;
  for(i=0; i<n; ++i)
  {
    for(j=0; j<m; ++j)
    {
      printf("%d ",arr[i][j]);
    }
    printf("\n");
  }
}

void RDL_print2DCharArray(char **arr, unsigned n, unsigned m)
{
  unsigned i,j;
  for(i=0; i<n; ++i)
  {
    for(j=0; j<m; ++j)
    {
      printf("%d ",arr[i][j]);
    }
    printf("\n");
  }
}

void RDL_writeNothing(RDL_ERROR_LEVEL level, const char* fmt, ...)
{

}

void RDL_writeToStderr(RDL_ERROR_LEVEL level, const char* fmt, ...)
{
  char* specifier;
  va_list va;

  /* we don't need any output initialization, so skip */
  if (level == RDL_INITIALIZE_OUTPUT) {
    return;
  }

  switch(level) {
    case RDL_DEBUG:
      specifier = "RDL_DEBUG";
      break;
    case RDL_WARNING:
      specifier = "RDL_WARNING";
      break;
    case RDL_ERROR:
      specifier = "RDL_ERROR";
      break;
    default:
      specifier = "RDL_????";
  }

  fprintf(stderr, "%12s: ", specifier);

  va_start(va, fmt);
  vfprintf(stderr, fmt, va);
  va_end(va);

}
