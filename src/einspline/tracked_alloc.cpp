
#include "tracked_alloc.h"
#include <Utilities/MemoryTracker_C.h>

void *maybe_align_alloc(size_t size, size_t alignment, const char *name)
{
  void *ptr;
#ifndef HAVE_POSIX_MEMALIGN
  ptr = malloc (size);
#else
  int err = posix_memalign (&ptr, alignment, size);
#endif
  MemoryTracker_add(ptr, size, name);
  return ptr;
}

void tracked_free(void *p)
{
  free(p);
  MemoryTracker_remove(p);
}
