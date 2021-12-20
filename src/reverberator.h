#ifndef REVERBERATOR_H
#define REVERBERATOR_H
#include "mainwindow.h"

class Reverberator
{
public:
    Reverberator();
    QVector<double> reverb(QVector<double> samples, float delay, float decayFactor, int mixPercent);
    QVector<double> combFilter(QVector<double> samples, int samplesLength, float delay, float decay, float sampleRate);
    QVector<double> allPassFilter(QVector<double> samples, int samplesLength, float sampleRate);
};

#endif // REVERBERATOR_H
