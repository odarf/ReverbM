#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "modeling.h"
#include "analysis.h"
#include "inou.h"
#include "qmath.h"

#include <iostream>
#include <QLogValueAxis>
#include <QLineSeries>
#include <QValueAxis>
#include <QChart>
#include <QChartView>
#include <QtCharts>
#include <QtWidgets>
#include <QtWidgets/QHBoxLayout>
#include <QHBoxLayout>
#include <cstdlib>

#include <math.h>
#include <time.h>
#include <random>

#include <fstream>

using namespace QtCharts;
using namespace std;

const int LENGTH = 1000;
const double DELTA = 50.0;

MainWindow::MainWindow(QWidget *parent): QMainWindow(parent), ui(new Ui::MainWindow) {
    ui->setupUi(this);
}

MainWindow::~MainWindow() { delete ui; }

void MainWindow::on_pushButton_clicked(){
    modeling  modeling;
    processing processing;
    analysis analysis;

    double k = -1;

    //double beta = 11;
    //double alpha = -0.01;

    int amplitude[] = {10, 100, 15};
    int frequency[] = {4, 37, 173};
    float randomCoeff = (float)ui->sbCoeff->value();

    QVector<double> x(LENGTH),             xLarge(3*LENGTH),
                    yCardio,
                    yMyRandom(LENGTH),     yEmbedRandom1(LENGTH),   yEmbedRandom2(LENGTH);

    yCardio.append(modeling.cardiogram());

    for(int X=0; X<LENGTH; X++){
        x[X] = X;
        yMyRandom[X] = modeling.randomGenerator(randomCoeff);
        yEmbedRandom1[X] = modeling.embedRandom();
        yEmbedRandom2[X] = modeling.embedRandom();
    }

    for(int i = 0; i<3000; i++){
        xLarge[i] = i;
    }

    double fc = 0.3;
    double fc2 = 0.6;
    int m = 64;
    int Loper = 2 * m + 1;
    QVector<double> lpfWeights = analysis.lowpassFilterPotter(fc, m);
    QVector<double> lpfWeights2 = analysis.lowpassFilterPotter(fc2, m);
    QVector<double> hpfWeights;
    QVector<double> pfWeights;
    QVector<double> rfWeights;

// ---------------- Рисование ----------------
    int currentTask = ui->tabWidget->currentIndex();
    ui->labelCI->setText(QString::number(currentTask));
    switch(currentTask){
        case 0:{
            QVector<double> dataFile = inou().load("C:/data.dat");
            QVector<double> xfile(dataFile.length());
            for (int i=0; i<xfile.length(); i++) { xfile[i] = i; }


            ui->graphFromFile->addGraph();
            ui->graphFromFile->graph(0)->setData(xfile, dataFile);
            ui->graphFromFile->xAxis->setRange(0, xfile.length());
            ui->graphFromFile->yAxis->setRange(analysis.minValue(dataFile)-10, analysis.maxValue(dataFile)+10);
            ui->graphFromFile->replot();

            QVector<double> fouramp = analysis.fourierAmplitude(dataFile);
            QVector<double> fourierFreq = analysis.calculateFrequency(0.002, dataFile.length());


            ui->graphFourierAmplitude->addGraph();
            ui->graphFourierAmplitude->graph(0)->setData(fourierFreq, fouramp);
            ui->graphFourierAmplitude->xAxis->setRange(0, 250);
            ui->graphFourierAmplitude->yAxis->setRange(analysis.minValue(fouramp), analysis.maxValue(fouramp)+5);
            ui->graphFourierAmplitude->replot();

            QVector<double> fourspec = analysis.fourierSpectrum(dataFile, 0.91);

            ui->graphFourierSpectrum->addGraph();
            ui->graphFourierSpectrum->graph(0)->setData(fourierFreq, fourspec);
            ui->graphFourierSpectrum->xAxis->setRange(0, 250);
            ui->graphFourierSpectrum->yAxis->setRange(analysis.minValue(fourspec), analysis.maxValue(fourspec)+6);
            ui->graphFourierSpectrum->replot();
            break;
        }
        case 1:
            ui->graphCardio->addGraph();
            ui->graphCardio->graph(0)->setData(x, yCardio);
            ui->graphCardio->xAxis->setRange(0, x.length());
            ui->graphCardio->yAxis->setRange(analysis.minValue(yCardio)-1, analysis.maxValue(yCardio)+1);
            ui->graphCardio->replot();
            break;
        case 2:{
            ui->graphLPF->addGraph();
            ui->graphLPF->graph(0)->setData(x, lpfWeights);
            ui->graphLPF->xAxis->setRange(0, 2*m);
            ui->graphLPF->yAxis->setRange(analysis.minValue(lpfWeights)-0.01, analysis.maxValue(lpfWeights)+0.01);
            ui->graphLPF->replot();

            for(int i = 0; i<=Loper-1; i++){ i==m ? hpfWeights.append(1 - lpfWeights[i]) : hpfWeights.append(-lpfWeights[i]); }

            ui->graphHPF->addGraph();
            ui->graphHPF->graph(0)->setData(x, hpfWeights);
            ui->graphHPF->xAxis->setRange(0, 2*m);
            ui->graphHPF->yAxis->setRange(analysis.minValue(hpfWeights)-0.01, analysis.maxValue(hpfWeights)+0.1);
            ui->graphHPF->replot();

            for(int i = 0; i<=Loper-1; i++){ pfWeights.append(lpfWeights2[i] - lpfWeights[i]); }

            ui->graphPF->addGraph();
            ui->graphPF->graph(0)->setData(x, pfWeights);
            ui->graphPF->xAxis->setRange(0, 2*m);
            ui->graphPF->yAxis->setRange(analysis.minValue(pfWeights)-0.01, analysis.maxValue(pfWeights)+0.01);
            ui->graphPF->replot();

            for(int i = 0; i<=Loper-1; i++){ i==m ? rfWeights.append(1 + lpfWeights[i] - lpfWeights2[i]) : rfWeights.append(lpfWeights[i] - lpfWeights2[i]); }

            ui->graphRF->addGraph();
            ui->graphRF->graph(0)->setData(x, rfWeights);
            ui->graphRF->xAxis->setRange(0, 2*m);
            ui->graphRF->yAxis->setRange(analysis.minValue(rfWeights)-0.01, analysis.maxValue(rfWeights)+0.01);
            ui->graphRF->replot();

            QVector<double> achxLP = analysis.fourierAmplitude(lpfWeights);
            QVector<double> achxHP = analysis.fourierAmplitude(hpfWeights);
            QVector<double> achxPF = analysis.fourierAmplitude(pfWeights);
            QVector<double> achxRF = analysis.fourierAmplitude(rfWeights);

            for(int i = 0; i<m; i++){
                achxLP[i] *= 2*m;
                achxHP[i] *= 2*m;
                achxPF[i] *= 2*m;
                achxRF[i] *= 2*m;
            }

            ui->graphAfcLp->addGraph();
            ui->graphAfcLp->graph(0)->setData(x, achxLP);
            ui->graphAfcLp->xAxis->setRange(0, m);
            ui->graphAfcLp->yAxis->setRange(0, analysis.maxValue(achxLP)+analysis.maxValue(achxLP)/10);
            ui->graphAfcLp->replot();

            ui->graphAfcHp->addGraph();
            ui->graphAfcHp->graph(0)->setData(x, achxHP);
            ui->graphAfcHp->xAxis->setRange(0, m);
            ui->graphAfcHp->yAxis->setRange(0, analysis.maxValue(achxHP)+analysis.maxValue(achxHP)/10);
            ui->graphAfcHp->replot();

            ui->graphAfcPf->addGraph();
            ui->graphAfcPf->graph(0)->setData(x, achxPF);
            ui->graphAfcPf->xAxis->setRange(0, m);
            ui->graphAfcPf->yAxis->setRange(0, analysis.maxValue(achxPF)+analysis.maxValue(achxPF)/10);
            ui->graphAfcPf->replot();

            ui->graphAfcRf->addGraph();
            ui->graphAfcRf->graph(0)->setData(x, achxRF);
            ui->graphAfcRf->xAxis->setRange(0, m);
            ui->graphAfcRf->yAxis->setRange(0, analysis.maxValue(achxRF)+analysis.maxValue(achxRF)/10);
            ui->graphAfcRf->replot();
            break;
        }
        case 3:{
            QVector<double> dataFile = inou().load("C:/data.dat");
            int N = dataFile.length();
            int M = lpfWeights.length();
            QVector<double> filteredLp(N+M);
            int tempN = N + M-1;

            for(auto i(0); i<tempN; ++i){
                int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
                int const jmx = (i < N - 1)  ? i           : N-1;
                for(auto j(jmn); j<=jmx; ++j){
                    filteredLp[i] += dataFile[j] * lpfWeights[i-j];
                }
            }

            for(int i(0); i<=M/2; i++){
                filteredLp.pop_front();
                filteredLp.pop_back();
            }

            ui->graphFilteredLp->addGraph();
            ui->graphFilteredLp->graph(0)->setData(x, filteredLp);
            ui->graphFilteredLp->xAxis->setRange(0, 1000);
            ui->graphFilteredLp->yAxis->setRange(analysis.minValue(filteredLp), analysis.maxValue(filteredLp)+analysis.maxValue(filteredLp)/10);
            ui->graphFilteredLp->replot();

            ui->graphFilteredLpAmpl->addGraph();
            ui->graphFilteredLpAmpl->graph(0)->setData(x, analysis.fourierAmplitude(filteredLp));
            ui->graphFilteredLpAmpl->xAxis->setRange(0, 400);
            ui->graphFilteredLpAmpl->yAxis->setRange(analysis.minValue(analysis.fourierAmplitude(filteredLp)), analysis.maxValue(analysis.fourierAmplitude(filteredLp))+analysis.maxValue(analysis.fourierAmplitude(filteredLp))/10);
            ui->graphFilteredLpAmpl->replot();

            M = hpfWeights.length();
            QVector<double> filteredHp(N+M);
            tempN = N + M-1;

            for(auto i(0); i<tempN; ++i){
                int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
                int const jmx = (i < N - 1)  ? i           : N-1;
                for(auto j(jmn); j<=jmx; ++j){
                    filteredHp[i] += dataFile[j] * hpfWeights[i-j];
                }
            }

            for(int i(0); i<=M/2; i++){
                filteredHp.pop_front();
                filteredHp.pop_back();
            }

            ui->graphFilteredHp->addGraph();
            ui->graphFilteredHp->graph(0)->setData(x, filteredHp);
            ui->graphFilteredHp->xAxis->setRange(0, 1000);
            ui->graphFilteredHp->yAxis->setRange(analysis.minValue(filteredHp), analysis.maxValue(filteredHp)+analysis.maxValue(filteredHp)/10);
            ui->graphFilteredHp->replot();

            ui->graphFilteredHpAmpl->addGraph();
            ui->graphFilteredHpAmpl->graph(0)->setData(x, analysis.fourierAmplitude(filteredHp));
            ui->graphFilteredHpAmpl->xAxis->setRange(0, 400);
            ui->graphFilteredHpAmpl->yAxis->setRange(analysis.minValue(analysis.fourierAmplitude(filteredHp)), analysis.maxValue(analysis.fourierAmplitude(filteredHp))+analysis.maxValue(analysis.fourierAmplitude(filteredHp))/10);
            ui->graphFilteredHpAmpl->replot();

            M = rfWeights.length();
            QVector<double> filteredRf(N+M);
            tempN = N + M-1;

            for(auto i(0); i<tempN; ++i){
                int const jmn = (i >= M - 1) ? i - (M - 1) : 0;
                int const jmx = (i < N - 1)  ? i           : N-1;
                for(auto j(jmn); j<=jmx; ++j){
                    filteredRf[i] += dataFile[j] * rfWeights[i-j];
                }
            }

            for(int i(0); i<=M/2; i++){
                filteredRf.pop_front();
                filteredRf.pop_back();
            }

            ui->graphFilteredPf->addGraph();
            ui->graphFilteredPf->graph(0)->setData(x, filteredRf);
            ui->graphFilteredPf->xAxis->setRange(0, 1000);
            ui->graphFilteredPf->yAxis->setRange(analysis.minValue(filteredRf), analysis.maxValue(filteredRf)+analysis.maxValue(filteredRf)/10);
            ui->graphFilteredPf->replot();

            ui->graphFilteredPfAmpl->addGraph();
            ui->graphFilteredPfAmpl->graph(0)->setData(x, analysis.fourierAmplitude(filteredRf));
            ui->graphFilteredPfAmpl->xAxis->setRange(0, 400);
            ui->graphFilteredPfAmpl->yAxis->setRange(analysis.minValue(analysis.fourierAmplitude(filteredRf)), analysis.maxValue(analysis.fourierAmplitude(filteredRf))+analysis.maxValue(analysis.fourierAmplitude(filteredRf))/10);
            ui->graphFilteredPfAmpl->replot();
            break;
        }
        case 4:{
        // --------------- Task 12 ----------------
            QVector<double> wav = inou().loadWave("../soundInput.wav");
            QVector<double> wavFourier = analysis.fourierAmplitude(wav);
            QVector<double> xWav(wav.length());
            auto min = analysis.minValue(wav);
            auto max = analysis.maxValue(wav);

            for (int i(0); i<wav.length(); i++) { xWav[i] = i; }

            ui->graphWavFile->addGraph();
            ui->graphWavFile->graph(0)->setData(xWav, wav);
            ui->graphWavFile->xAxis->setRange(0, xWav.length());
            ui->graphWavFile->yAxis->setRange(min, max+max/10);
            ui->graphWavFile->replot();

            min = analysis.minValue(wavFourier);
            max = analysis.maxValue(wavFourier);

            ui->graphWavFileFourier->addGraph();
            ui->graphWavFileFourier->graph(0)->setData(xWav, wavFourier);
            ui->graphWavFileFourier->xAxis->setRange(0, wavFourier.length());
            ui->graphWavFileFourier->yAxis->setRange(min, max+max/10);
            ui->graphWavFileFourier->replot();

            //for (int i(0); i<wav.length(); i++) { wav[i] *= 4; }

            min = analysis.minValue(wav);
            max = analysis.maxValue(wav);

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(min, max+max/10);
            ui->graphWavFileLoud->replot();
        // ----------------------------------------
            break;
        }
        case 5:{
            //ПРОДОЛЖИТЬ
            QVector<double> wav = inou().loadWave("../soundInput.wav");
            QVector<double> Fonem;

            for(int i(4000); i<16000; i++){
                Fonem.append(wav[i]);
            }
            QVector<double> xWav;
            for (int i(0); i<wav.length(); i++) {
                xWav.append(i);
            }


            QVector<double> FonemFourier = analysis.fourierAmplitude(Fonem);

            ui->graph08->addGraph();
            ui->graph08->graph(0)->setData(xWav, Fonem);
            ui->graph08->xAxis->setRange(0, Fonem.length());
            ui->graph08->yAxis->setRange(analysis.minValue(Fonem)-analysis.maxValue(Fonem)/10, analysis.maxValue(Fonem)+analysis.maxValue(Fonem)/10);
            ui->graph08->replot();

            int N = Fonem.length();
            //double dt = 1 / FonemFourier.length(); //
            //if(dt = 0){ throw std::exception(); }
            double f_gr = N/2;
            double df = (2 * f_gr) / N;
            //float df = 22050/FonemFourier.length();

            QVector<double> xWav1;
            for (int i= 0; i<FonemFourier.length(); i++) {
                xWav1.append((double)(i * df));
            }
            inou().exportWave(Fonem, Fonem.length(), "../123.wav", 1);

            ui->graph08_2->addGraph();
            ui->graph08_2->graph(0)->setData(xWav1, FonemFourier);
            ui->graph08_2->xAxis->setRange(0, FonemFourier.length());
            ui->graph08_2->setInteraction(QCP::iRangeZoom,true);
            ui->graph08_2->setInteraction(QCP::iRangeDrag,true);
            ui->graph08_2->yAxis->setRange(analysis.minValue(FonemFourier)-analysis.maxValue(FonemFourier)/10, analysis.maxValue(FonemFourier)+analysis.maxValue(FonemFourier)/10);
            ui->graph08_2->replot();

            /*ui->graph08_3->addGraph();
            ui->graph08_3->graph(0)->setData(xWav, foo);
            ui->graph08_3->xAxis->setRange(0, foo.length());
            ui->graph08_3->yAxis->setRange(analysis.minValue(foo)-analysis.maxValue(foo)/10, analysis.maxValue(foo)+analysis.maxValue(foo)/10);
            ui->graph08_3->replot();*/

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(analysis.minValue(wav)-analysis.maxValue(wav)/10, analysis.maxValue(wav)+analysis.maxValue(wav)/10);
            ui->graphWavFileLoud->replot();

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(analysis.minValue(wav)-analysis.maxValue(wav)/10, analysis.maxValue(wav)+analysis.maxValue(wav)/10);
            ui->graphWavFileLoud->replot();

            ui->graphWavFileLoud->addGraph();
            ui->graphWavFileLoud->graph(0)->setData(xWav, wav);
            ui->graphWavFileLoud->xAxis->setRange(0, xWav.length());
            ui->graphWavFileLoud->yAxis->setRange(analysis.minValue(wav)-analysis.maxValue(wav)/10, analysis.maxValue(wav)+analysis.maxValue(wav)/10);
            ui->graphWavFileLoud->replot();
        }
        default:
            return;
    }

// ---------- Очистка графиков ----------

// --------------------------------------
}
