#include <QCoreApplication>

#include "gifgenerator.h"

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);
    GifGenerator gifgenerator(QString("./"), 600, 800, 100);
    QImage img1("./1.png");
    QImage img2("./2.png");
    QImage img3("./3.png");

    gifgenerator.add(img1);
    gifgenerator.add(img2);
    gifgenerator.add(img2);
    gifgenerator.save();
    return a.exec();
}
