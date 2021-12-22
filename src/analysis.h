#ifndef ANALYSIS_H
#define ANALYSIS_H

#include <QMainWindow>

class analysis
{
private:
public:
    analysis();

    double minValue(QVector<double> x);

    double maxValue(QVector<double> x);

    QVector<double> fourierAmplitude(QVector<double> inputData);


};

#endif // ANALYSIS_H
