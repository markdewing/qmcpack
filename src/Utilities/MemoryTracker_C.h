
#ifndef MEMORYTRACKER_C
#define MEMORYTRACKER_C

#ifdef __cplusplus
extern "C" {
#endif
  void MemoryTracker_add(void *ptr, size_t size, const char *name);
  void MemoryTracker_remove(void *ptr);
#ifdef __cplusplus
}
#endif

#endif
