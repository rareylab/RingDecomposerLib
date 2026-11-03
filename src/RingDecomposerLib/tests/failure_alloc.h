/* Allocator hooks used only by the fault-injection build. */
#ifndef RDL_FAILURE_ALLOC_H
#define RDL_FAILURE_ALLOC_H
#include <stddef.h>
void *audit_malloc(size_t);
void *audit_calloc(size_t,size_t);
void *audit_realloc(void *,size_t);
void audit_free(void *);
#endif
