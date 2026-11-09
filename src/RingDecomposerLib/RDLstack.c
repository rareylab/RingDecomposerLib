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

#include "RDLstack.h"

#include <stdlib.h>

#include "RDLutility.h"
#include <limits.h>

struct RDL_stack {
  void **elements;
  unsigned size, reserved;
};

RDL_stack *RDL_stack_new(void)
{
  RDL_stack *stack = RDL_alloc(1, sizeof(*stack));
  if (!stack) {
    return NULL;
  }
  stack->reserved = 32;
  stack->elements = RDL_alloc(stack->reserved, sizeof(*stack->elements));
  if (!stack->elements) { RDL_free(stack); return NULL; }
  return stack;
}

void *RDL_stack_top(RDL_stack *stack)
{
  return stack->size ? stack->elements[stack->size - 1] : NULL;
}

int RDL_stack_empty(RDL_stack *stack) { return stack->size == 0; }
void RDL_stack_pop(RDL_stack *stack)
{
  if (stack->size) {
    --stack->size;
  }
}

int RDL_stack_push(RDL_stack *stack, void *element)
{
  void **elements;
  if (stack->size == stack->reserved) {
    if (stack->reserved > UINT_MAX / 2) {
      return 0;
    }
    elements = RDL_resize(stack->elements, (size_t)stack->reserved * 2,
                          sizeof(*elements));
    if (!elements) {
      return 0;
    }
    stack->elements = elements;
    stack->reserved *= 2;
  }
  stack->elements[stack->size++] = element;
  return 1;
}

void RDL_stack_delete(RDL_stack *stack)
{
  if (!stack) {
    return;
  }
  RDL_free(stack->elements);
  RDL_free(stack);
}
