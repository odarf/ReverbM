#include "reverberator.h"

Reverberator::Reverberator(){ }

QVector<double> Reverberator::reverb(QVector<double> samples, float delay, float decayFactor, int mixPercent){
    const float sampleRate = 22.050f;
    int length = samples.length();

    QVector<double> combFilterSamples1 = combFilter(samples, length, delay, decayFactor, sampleRate);
    QVector<double> combFilterSamples2 = combFilter(samples, length, delay - 11.73f, decayFactor, sampleRate);
    QVector<double> combFilterSamples3 = combFilter(samples, length, delay + 19.31f, decayFactor, sampleRate);
    QVector<double> combFilterSamples4 = combFilter(samples, length, delay - 7.97f, decayFactor, sampleRate);

    QVector<double> outputComb(length);

    for(int i = 0; i<length; i++){
        outputComb[i] = combFilterSamples1[i] + combFilterSamples2[i] + combFilterSamples3[i] + combFilterSamples4[i];
    }

    combFilterSamples1.clear();
    combFilterSamples2.clear();
    combFilterSamples3.clear();
    combFilterSamples4.clear();


    QVector<double> mix(length);
    for(int i = 0; i<length; i++){
        mix[i] = ((100 - mixPercent) * samples[i]) + (mixPercent * outputComb[i]);
    }

    QVector<double> allPassFilterSamples1 = allPassFilter(mix, length, sampleRate);
    QVector<double> allPassFilterSamples2 = allPassFilter(allPassFilterSamples1, length, sampleRate);

    return allPassFilterSamples2;
}

QVector<double> Reverberator::combFilter(QVector<double> samples, int samplesLength, float delay, float decay, float sampleRate){
    int delaySamples = (int)((float)delay * 22.050f);
    QVector<double> combFilterSamples = samples;

    for(int i = 0; i<samplesLength - delaySamples; i++){
        combFilterSamples[i+delaySamples] +=((float)combFilterSamples[i] * decay);
    }
    return combFilterSamples;
}

QVector<double> Reverberator::allPassFilter(QVector<double> samples, int samplesLength, float sampleRate){
    //int delaySamples = (int)((float)89.27f * (sampleRate/1000));
    int delaySamples = (int)((float)89.27f * 22.050f);
    QVector<double> allPassFilterSamples(samplesLength);
    float decayFactor = 0.131f;

    for(int i = 0; i<samplesLength; i++){
        allPassFilterSamples[i] = samples[i];

        if(i-delaySamples >= 0){
            allPassFilterSamples[i] += -decayFactor * allPassFilterSamples[i-delaySamples];
        }

        if(i-delaySamples >= 1){
            allPassFilterSamples[i] += decayFactor * allPassFilterSamples[i+20-delaySamples];
        }
    }
    //Нормализация до 1
    float value = allPassFilterSamples[0];
    float max = 0.0f;

    for(int i = 0; i<samplesLength; i++){
        if(std::abs(allPassFilterSamples[i]) > max){
            max = std::abs(allPassFilterSamples[i]);
        }
    }

    for(int i = 0; i<allPassFilterSamples.length(); i++){
        float currentValue = allPassFilterSamples[i];
        value = ((value + (currentValue - value)) / max);
        allPassFilterSamples[i] = value;
    }

    return allPassFilterSamples;
}
