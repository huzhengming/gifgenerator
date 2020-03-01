#ifndef GIFGENERATOR_H
#define GIFGENERATOR_H

/*
作者：胡正明
联系方式：8615827082382
邮箱：1959899835@qq.com

代码功能：能够根据不同的图片生成gif 动图

主要方便QT开发者理解gif 的编码原理，和简单实用

使用示例
  GifGenerator gifgenerator(QString("./"), 600, 800);
  QImage img("./1.png");
  QImage img2("./2.png");

  gifgenerator.add(img);
  gifgenerator.add(img2);
  gifgenerator.save();
*/

#include <QObject>
#include <QString>
#include <QImage>

#include <stdlib.h>


class GifGenerator : public QObject
{
    Q_OBJECT
public:
    struct Palette
    {
        int bitDepth;

        uint8_t r[256];
        uint8_t g[256];
        uint8_t b[256];

        // k-d tree over RGB space, organized in heap fashion
        // i.e. left child of node i is node i*2, right child is node i*2+1
        // nodes 256-511 are implicitly the leaves, containing a color
        uint8_t treeSplitElt[255];
        uint8_t treeSplit[255];
    };

    // The LZW dictionary is a 256-ary tree constructed as the file is encoded,
    // this is one node
    struct LzwNode
    {
        uint16_t m_next[256];
    };

    // Simple structure to write out the LZW-compressed portion of the image
    // one bit at a time
    struct BitStatus
    {
        uint8_t bitIndex;  // how many bits in the partial byte written so far
        uint8_t byte;      // current partial byte

        int chunkIndex;
        uint8_t chunk[256];   // bytes are written in here until we have 256 of them, then written to the file
    };



    explicit GifGenerator(QString gifPath, int width, int height, int delay = 40 ,QObject *parent = nullptr);
    ~ GifGenerator();
    bool init(QString gifPath);
    bool reset(QString gifPath, int width, int height);
    bool save();
    bool add(QImage & img);
protected:
    void MakePalette( const uint8_t* lastFrame, const uint8_t* nextFrame, int bitDepth, bool buildForDither, Palette* pPal );
    void DitherImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, Palette* pPal );
    void ThresholdImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, Palette* pPal );
    void WriteLzwImage(uint8_t* image, int left, int top, int delay, Palette* pPal);
    int  PickChangedPixels( const uint8_t* lastFrame, uint8_t* frame, int numPixels );
    void SplitPalette(uint8_t* image, int numPixels, int firstElt, int lastElt, int splitElt, int splitDist, int treeNode, bool buildForDither, Palette* pal);
    void GetClosestPaletteColor(Palette* pPal, int r, int g, int b, int& bestInd, int& bestDiff, int treeRoot = 1);
    void SwapPixels(uint8_t* image, int pixA, int pixB);
    int Partition(uint8_t* image, const int left, const int right, const int elt, int pivotIndex);
    void PartitionByMedian(uint8_t* image, int left, int right, int com, int neededCenter);
    void WriteBit( BitStatus& stat, int bit );
    void WriteChunk(BitStatus& stat );
    void WriteCode(BitStatus& stat, int code, int length );
    void WritePalette( const Palette* pPal);
    // max, min, and abs functions
    int IMax(int l, int r) { return l>r?l:r; }
    int IMin(int l, int r) { return l<r?l:r; }
    int IAbs(int i) { return i<0?-i:i; }
public slots:
private:
    bool m_initOk;
    int m_width;
    int m_height;
    int m_delay;  //帧间隔 40 ms
    bool m_firstFrame;
    int m_kGifTransIndex;
    QString m_gifPath;
    FILE* m_fpGif;
    uint8_t* m_oldImage;


};

#endif // GIFGENERATOR_H
