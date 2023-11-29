#include <stdio.h>
#include <thread>
#include <iostream>
#include <vector>

using std::vector;
using std::cout;
using std::endl;

#include "CycleTimer.h"

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
    int startRow;
    int numRows;
    double time;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);

// re: 实部; im: 虚部
// 判断某点是否属于mandel集合, 返回迭代次数
static inline int mandel(float c_re, float c_im, int count)
{
    float z_re = c_re, z_im = c_im;
    int i;
    for (i = 0; i < count; ++i) {
        // 判断Mandelbrot集合的迭代过程中，当前的复数是否已经趋于无穷大。
        // Mandelbrot集合中的点在迭代过程中要么趋于无穷大，要么保持有界（复平面上，Mandelbrot集合的点在迭代过程中始终保持在一个半径为2的圆内）
        if (z_re * z_re + z_im * z_im > 4.f)
            break;
        // iteration
        float new_re = z_re*z_re - z_im*z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }

    return i;
}

//
// workerThreadStart --
//
// Thread entrypoint.
void workerThreadStart(WorkerArgs * const args) {

    // TODO FOR CS149 STUDENTS: Implement the body of the worker
    // thread here. Each thread should make a call to mandelbrotSerial()
    // to compute a part of the output image.  For example, in a
    // program that uses two threads, thread 0 could compute the top
    // half of the image and thread 1 could compute the bottom half.

    printf("Hello world from thread %d\n", args->threadId);

    double startTime = CycleTimer::currentSeconds();

    float dx = (args->x1 - args->x0) / args->width;                 // 每个小格子的宽
    float dy = (args->y1 - args->y0) / args->height;                // 每个小格子的高
    int endRow = args->startRow + args->numRows;
    for (int j = args->startRow; j < endRow; ++j) {
        for (unsigned int i = 0; i < args->width; ++i) {
            float x = args->x0 + i * dx;
            float y = args->y0 + j * dy;
            int index = j * args->width + i ;
            args->output[index] = mandel(x, y, args->maxIterations);
        }
    }

    double endTime = CycleTimer::currentSeconds();
    args->time = (endTime - startTime) * 1000;
}


//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
vector<double> mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    static constexpr int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    vector<double> vec(numThreads);

    // Creates thread objects that do not yet represent a thread.
    std::thread workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];
    
    int y_interval = height / numThreads;     //  每个小矩形的height有多少个格子
    for (int i=0; i<numThreads; i++) {
        // TODO FOR CS149 STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
        args[i].threadId = i;
        if (i == numThreads - 1) {                  // 本线程负责多少行
            args[i].numRows = height - y_interval * i;
        } else {
            args[i].numRows = y_interval;
        }
        args[i].startRow = i * y_interval;          // 起始行
    }

    // Spawn the worker threads.  Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
    for (int i=1; i<numThreads; i++) {
        workers[i] = std::thread(workerThreadStart, &args[i]);
    }
    
    workerThreadStart(&args[0]);

    // join worker threads
    for (int i=1; i<numThreads; i++) {
        workers[i].join();
    }

    for (int i = 0; i < numThreads; ++i) {
        vec[i] = args[i].time;
    }

    return vec;
}

