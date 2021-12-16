#ifndef MODELING_H
#define MODELING_H

#include <iostream>
#include <QMainWindow>


class modeling
{
public:
    modeling();

    ///Встроенный рандом(?)
    float embedRandom();

    ///Мой рандом(?)
    float randomGenerator(float coefficient);

    ///Модель кардиограммы
    QVector<double> cardiogram();
};

#endif // MODELING_H
