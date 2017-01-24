
#ifndef TRACKED_ALLOC_H
#define TRACKED_ALLOC_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

int posix_memalign(void **memptr, size_t alignment, size_t size);

void *maybe_align_alloc(size_t size, size_t alignment, const char *name);
#if 0
{
  void *ptr;
#ifndef HAVE_POSIX_MEMALIGN
  ptr = malloc (sizeof(float)*Nx*N);
#else
  int err = posix_memalign (&ptr, alignment, size);
#endif
  MemoryTracker_add(ptr, size, name);
  return ptr;
}
#endif

void tracked_free(void *p);

#ifdef __cplusplus
}
#endif
#endif

