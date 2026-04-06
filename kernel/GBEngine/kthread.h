// Stub header for parallel bba threading support
// Not yet implemented — this file exists only to satisfy #include
#ifndef KTHREAD_H
#define KTHREAD_H

#include "kernel/GBEngine/kutil.h"

struct SweepContext;

static inline int get_singular_threads() {
  const char *env = getenv("SINGULAR_THREADS");
  if (env != NULL) return atoi(env);
  return 1;
}

static inline SweepContext* sweep_context_init(kStrategy /*strat*/, int /*threads*/) {
  return NULL;
}

static inline void bba_parallel_loop(SweepContext* /*ctx*/) {
}

static inline void sweep_context_destroy(SweepContext* /*ctx*/) {
}

#endif
