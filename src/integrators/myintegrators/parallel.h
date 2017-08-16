#pragma once

#include <mitsuba/mitsuba.h>
#include <mutex>
#include <functional>
#include <atomic>
#include <vector>

// From https://github.com/mmp/pbrt-v3/blob/master/src/core/parallel.h

MTS_NAMESPACE_BEGIN

#if defined(__APPLE__) && defined(__MACH__)
// clang on MacOS does not support thread_local keyword...
#define thread_local __thread
#endif

void ParallelFor(const std::function<void(int64_t)> &func, int64_t count, int chunkSize = 1);
extern thread_local int threadIndex;
void ParallelFor(std::function<void(Vector2i)> func, const Vector2i count);
int MaxThreadIndex();
int NumSystemCores();
void TerminateWorkerThreads();

MTS_NAMESPACE_END
