#include "modeling.h"
#include "qmath.h"

#include <iostream>
#include <QMainWindow>

#include <cstdlib>

#include <math.h>
#include <random>
#include <time.h>

using namespace std;

modeling::modeling() {}

float modeling::embedRandom(){
    int const max = 1;
    int const min = -1;
    return static_cast<float>(rand() * (1.0 / (static_cast<float>(RAND_MAX) + 1.0)) * (max - min) + min);
}

float modeling::randomGenerator(float coefficient){
    float const left = -1.0 * coefficient;
    float const right = 1.0 * coefficient;
    default_random_engine engine{random_device()()};
    uniform_real_distribution<float> distribution{left, right};
    return distribution(engine);
}

QVector<double> modeling::cardiogram(){
    int const N = 1000;
    int const M = 200;
    QVector<double> x(N), h(M), y(N+M);

    float alpha = 10;
    int frequency = 4;
    int const tempN = N + M - 1;

    for(int i = 0; i<h.length(); i++){
        h[i] = sin(2 * 3.14 * frequency * (i*0.005)) * exp(-alpha * (i*0.005));
    }

    x[200] = 120;
    x[400] = 130;
    x[600] = 110;

    for(auto i(0); i<tempN; ++i){
        int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
        int const jmx = (i < N - 1)  ? i           : N-1;
        for(auto j(jmn); j<=jmx; ++j){
            y[i] += x[j] * h[i-j];
        }
    }
    return y;
}
