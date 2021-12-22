#include "analysis.h"
#include "mainwindow.h"

#include <algorithm>
#include <math.h>
#include <QDebug>

analysis::analysis(){}

double analysis::minValue(QVector<double> x){
    double min = *std::min_element(x.constBegin(), x.constEnd());
    return min;
}

double analysis::maxValue(QVector<double> x){
    double max = *std::max_element(x.constBegin(), x.constEnd());
    return max;
}

QVector<double> analysis::fourierAmplitude(QVector<double> inputData){
    int length = inputData.length();
    QVector<double> outputData;

    for(int i = 0; i<length/2; i++){
        double real = 0;
        double imagine = 0;

        for(int j = 0; j<length; j++){
            real += inputData[j] * std::cos((2 * 3.14 * i * j) / length);
            imagine += inputData[j] * std::sin((2 * 3.14 * i * j) / length);
        }
        real /= length;
        imagine /= length;
        outputData.append(std::sqrt(real*real + imagine*imagine));
    }
    return outputData;
}
