
#ifndef GLIN_STOPWATCH_H
#define GLIN_STOPWATCH_H
#include <chrono>

class Stopwatch {
private:
    std::chrono::high_resolution_clock::time_point t1, t2;

public:
    explicit Stopwatch(bool run = false) {
        if (run) {
            start();
        }
    }

    void start() { t2 = t1 = std::chrono::high_resolution_clock::now(); }
    void stop() { t2 = std::chrono::high_resolution_clock::now(); }

    double ms() const {
        return std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1)
                       .count() /
               1000.0;
    }
};
#endif //GLIN_STOPWATCH_H
