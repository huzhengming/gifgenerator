#include "gifgenerator.h"

#include <QDateTime>
#include <QDir>

GifGenerator::GifGenerator(QString gifPath, int width, int height, int delay ,QObject * parent )
    : QObject(parent)
    , m_width(width)
    , m_height(height)
    , m_kGifTransIndex(0)
    , m_oldImage(nullptr)
    , m_firstFrame(true)
    , m_delay(delay)
{
    init(gifPath);
}

GifGenerator:: ~ GifGenerator()
{
    if (m_initOk && m_fpGif) {
       fclose(m_fpGif);
    }
}

bool GifGenerator::init(QString gifPath)
{
    m_initOk = true;
    //判断文件路径是否合法
    m_fpGif = fopen(gifPath.toStdString().c_str(),"wb+");
    if (!m_fpGif) {
        QDir dir(gifPath);
        QString pathHead = "";
        if (dir.exists()) {
            if(gifPath.back() != '/') {
                gifPath += '/';
            }
            pathHead = gifPath;
        }

        QDateTime time = QDateTime::currentDateTime();
        QString path = pathHead + time.toString("yyyy_MM_dd_hh_mm_ss") + ".gif";
        m_fpGif = fopen(path.toStdString().c_str(),"wb+");
        if (!m_fpGif) {
            m_initOk = false;
        }
        else {
            m_gifPath = path;
        }
    }
    else {
        m_gifPath = gifPath;
    }

    if (m_initOk) {
        // allocate
        m_oldImage = (uint8_t*)malloc(m_width*m_height*4);

        fputs("GIF89a", m_fpGif);

        // screen descriptor
        fputc(m_width & 0xff, m_fpGif);
        fputc((m_width >> 8) & 0xff, m_fpGif);
        fputc(m_height & 0xff, m_fpGif);
        fputc((m_height >> 8) & 0xff, m_fpGif);

        fputc(0xf0, m_fpGif);  // there is an unsorted global color table of 2 entries
        fputc(0, m_fpGif);     // background color
        fputc(0, m_fpGif);     // pixels are square (we need to specify this because it's 1989)

        // now the "global" palette (really just a dummy palette)
        // color 0: black
        fputc(0, m_fpGif);
        fputc(0, m_fpGif);
        fputc(0, m_fpGif);
        // color 1: also black
        fputc(0, m_fpGif);
        fputc(0, m_fpGif);
        fputc(0, m_fpGif);

        //========delay != 0
        // animation header
        fputc(0x21, m_fpGif); // extension
        fputc(0xff, m_fpGif); // application specific
        fputc(11, m_fpGif); // length 11
        fputs("NETSCAPE2.0", m_fpGif); // yes, really
        fputc(3, m_fpGif); // 3 bytes of NETSCAPE2.0 data

        fputc(1, m_fpGif); // JUST BECAUSE
        fputc(0, m_fpGif); // loop infinitely (byte 0)
        fputc(0, m_fpGif); // loop infinitely (byte 1)

        fputc(0, m_fpGif); // block terminator
    }
    return m_initOk;
}

bool GifGenerator::reset(QString gifPath, int width, int height)
{
    if (m_initOk && m_fpGif) {
       fclose(m_fpGif);
       m_fpGif = nullptr;
    }
     m_width = width;
     m_height = height;
    return init(gifPath);
}

bool GifGenerator::save()
{
    if (!m_initOk || !m_fpGif) return false;

    fputc(0x3b, m_fpGif); // end of file
    fclose(m_fpGif);
    free(m_oldImage);

    m_fpGif = nullptr;
    m_oldImage = nullptr;

    return true;
}

bool GifGenerator::add(QImage & img)
{
    if (!m_initOk || !m_fpGif) return false;
    if (img.width() == 0 || img.height() == 0) return false;
    QImage image= img.scaled(m_width,m_height).convertToFormat(QImage::Format_RGBA8888);
    const uint8_t* oldImage = m_firstFrame? nullptr : m_oldImage;
    m_firstFrame = false;

    bool dither = false;
    int bitDepth = 8;//像素深度 8   
    Palette pal;
    MakePalette((dither? nullptr : oldImage), image.bits(),bitDepth, dither, &pal);

    if (dither)
        DitherImage(oldImage, image.bits(), m_oldImage, &pal);
    else
        ThresholdImage(oldImage, image.bits(), m_oldImage, &pal);

    WriteLzwImage( m_oldImage, 0, 0, m_delay, &pal);

    return true;
}

void GifGenerator::MakePalette( const uint8_t* lastFrame, const uint8_t* nextFrame, int bitDepth, bool buildForDither, Palette* pPal )
{
    pPal->bitDepth = bitDepth;

    // SplitPalette is destructive (it sorts the pixels by color) so
    // we must create a copy of the image for it to destroy
    size_t imageSize = (size_t)(m_width * m_height * 4 * sizeof(uint8_t));
    uint8_t* destroyableImage = (uint8_t*)malloc(imageSize);
    memcpy(destroyableImage, nextFrame, imageSize);

    int numPixels = (int)(m_width * m_height);
    if(lastFrame)
        numPixels = PickChangedPixels(lastFrame, destroyableImage, numPixels);

    const int lastElt = 1 << bitDepth;
    const int splitElt = lastElt/2;
    const int splitDist = splitElt/2;

    SplitPalette(destroyableImage, numPixels, 1, lastElt, splitElt, splitDist, 1, buildForDither, pPal);

    free(destroyableImage);

    // add the bottom node for the transparency index
    pPal->treeSplit[1 << (bitDepth-1)] = 0;
    pPal->treeSplitElt[1 << (bitDepth-1)] = 0;

    pPal->r[0] = pPal->g[0] = pPal->b[0] = 0;
}

void GifGenerator::DitherImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame,  Palette* pPal )
{
    int numPixels = (int)(m_width * m_height);

    // quantPixels initially holds color*256 for all pixels
    // The extra 8 bits of precision allow for sub-single-color error values
    // to be propagated
    int32_t *quantPixels = (int32_t *)malloc(sizeof(int32_t) * (size_t)numPixels * 4);

    for ( int ii=0; ii<numPixels*4; ++ii )
    {
        uint8_t pix = nextFrame[ii];
        int32_t pix16 = int32_t(pix) * 256;
        quantPixels[ii] = pix16;
    }

    for ( uint32_t yy=0; yy<m_height; ++yy )
    {
        for( uint32_t xx=0; xx<m_width; ++xx )
        {
            int32_t* nextPix = quantPixels + 4*(yy*m_width+xx);
            const uint8_t* lastPix = lastFrame? lastFrame + 4*(yy*m_width+xx) : nullptr;

            // Compute the colors we want (rounding to nearest)
            int32_t rr = (nextPix[0] + 127) / 256;
            int32_t gg = (nextPix[1] + 127) / 256;
            int32_t bb = (nextPix[2] + 127) / 256;

            // if it happens that we want the color from last frame, then just write out
            // a transparent pixel
            if( lastFrame &&
               lastPix[0] == rr &&
               lastPix[1] == gg &&
               lastPix[2] == bb )
            {
                nextPix[0] = rr;
                nextPix[1] = gg;
                nextPix[2] = bb;
                nextPix[3] = m_kGifTransIndex;
                continue;
            }

            int32_t bestDiff = 1000000;
            int32_t bestInd = m_kGifTransIndex;

            // Search the palete
            GetClosestPaletteColor(pPal, rr, gg, bb, bestInd, bestDiff);

            // Write the result to the temp buffer
            int32_t r_err = nextPix[0] - int32_t(pPal->r[bestInd]) * 256;
            int32_t g_err = nextPix[1] - int32_t(pPal->g[bestInd]) * 256;
            int32_t b_err = nextPix[2] - int32_t(pPal->b[bestInd]) * 256;

            nextPix[0] = pPal->r[bestInd];
            nextPix[1] = pPal->g[bestInd];
            nextPix[2] = pPal->b[bestInd];
            nextPix[3] = bestInd;

            // Propagate the error to the four adjacent locations
            // that we haven't touched yet
            int quantloc_7 = (int)(yy * m_width + xx + 1);
            int quantloc_3 = (int)(yy * m_width + m_width + xx - 1);
            int quantloc_5 = (int)(yy * m_width + m_width + xx);
            int quantloc_1 = (int)(yy * m_width + m_width + xx + 1);

            if(quantloc_7 < numPixels)
            {
                int32_t* pix7 = quantPixels+4*quantloc_7;
                pix7[0] += IMax( -pix7[0], r_err * 7 / 16 );
                pix7[1] += IMax( -pix7[1], g_err * 7 / 16 );
                pix7[2] += IMax( -pix7[2], b_err * 7 / 16 );
            }

            if(quantloc_3 < numPixels)
            {
                int32_t* pix3 = quantPixels+4*quantloc_3;
                pix3[0] += IMax( -pix3[0], r_err * 3 / 16 );
                pix3[1] += IMax( -pix3[1], g_err * 3 / 16 );
                pix3[2] += IMax( -pix3[2], b_err * 3 / 16 );
            }

            if(quantloc_5 < numPixels)
            {
                int32_t* pix5 = quantPixels+4*quantloc_5;
                pix5[0] += IMax( -pix5[0], r_err * 5 / 16 );
                pix5[1] += IMax( -pix5[1], g_err * 5 / 16 );
                pix5[2] += IMax( -pix5[2], b_err * 5 / 16 );
            }

            if(quantloc_1 < numPixels)
            {
                int32_t* pix1 = quantPixels+4*quantloc_1;
                pix1[0] += IMax( -pix1[0], r_err / 16 );
                pix1[1] += IMax( -pix1[1], g_err / 16 );
                pix1[2] += IMax( -pix1[2], b_err / 16 );
            }
        }
    }

    // Copy the palettized result to the output buffer
    for( int ii=0; ii<numPixels*4; ++ii )
    {
        outFrame[ii] = (uint8_t)quantPixels[ii];
    }

    free(quantPixels);
}

void GifGenerator::ThresholdImage( const uint8_t* lastFrame, const uint8_t* nextFrame, uint8_t* outFrame, Palette* pPal )
{
    uint32_t numPixels = m_width * m_height;
    for ( uint32_t ii=0; ii<numPixels; ++ii )
    {
        // if a previous color is available, and it matches the current color,
        // set the pixel to transparent
        if (lastFrame &&
           lastFrame[0] == nextFrame[0] &&
           lastFrame[1] == nextFrame[1] &&
           lastFrame[2] == nextFrame[2])
        {
            outFrame[0] = lastFrame[0];
            outFrame[1] = lastFrame[1];
            outFrame[2] = lastFrame[2];
            outFrame[3] = m_kGifTransIndex;
        }
        else
        {
            // palettize the pixel
            int32_t bestDiff = 1000000;
            int32_t bestInd = 1;
            GetClosestPaletteColor(pPal, nextFrame[0], nextFrame[1], nextFrame[2], bestInd, bestDiff);

            // Write the resulting color to the output buffer
            outFrame[0] = pPal->r[bestInd];
            outFrame[1] = pPal->g[bestInd];
            outFrame[2] = pPal->b[bestInd];
            outFrame[3] = (uint8_t)bestInd;
        }

        if(lastFrame) lastFrame += 4;
        outFrame += 4;
        nextFrame += 4;
    }
}

void GifGenerator::WriteLzwImage(uint8_t* image, int left, int top, int delay, Palette* pPal)
{
    // graphics control extension
    fputc(0x21, m_fpGif);
    fputc(0xf9, m_fpGif);
    fputc(0x04, m_fpGif);
    fputc(0x05, m_fpGif); // leave prev frame in place, this frame has transparency
    fputc(delay & 0xff, m_fpGif);
    fputc((delay >> 8) & 0xff, m_fpGif);
    fputc(m_kGifTransIndex, m_fpGif); // transparent color index
    fputc(0, m_fpGif);

    fputc(0x2c, m_fpGif); // image descriptor block

    fputc(left & 0xff, m_fpGif);           // corner of image in canvas space
    fputc((left >> 8) & 0xff, m_fpGif);
    fputc(top & 0xff, m_fpGif);
    fputc((top >> 8) & 0xff, m_fpGif);

    fputc(m_width & 0xff, m_fpGif);          // width and height of image
    fputc((m_width >> 8) & 0xff, m_fpGif);
    fputc(m_height & 0xff, m_fpGif);
    fputc((m_height >> 8) & 0xff, m_fpGif);

    //fputc(0, f); // no local color table, no transparency
    //fputc(0x80, f); // no local color table, but transparency

    fputc(0x80 + pPal->bitDepth-1, m_fpGif); // local color table present, 2 ^ bitDepth entries
    WritePalette(pPal);

    const int minCodeSize = pPal->bitDepth;
    const uint32_t clearCode = 1 << pPal->bitDepth;

    fputc(minCodeSize, m_fpGif); // min code size 8 bits

    LzwNode* codetree = (LzwNode*)malloc(sizeof(LzwNode)*4096);

    memset(codetree, 0, sizeof(LzwNode)*4096);
    int32_t curCode = -1;
    uint32_t codeSize = (uint32_t)minCodeSize + 1;
    uint32_t maxCode = clearCode+1;

    BitStatus stat;
    stat.byte = 0;
    stat.bitIndex = 0;
    stat.chunkIndex = 0;

    WriteCode(stat, clearCode, codeSize);  // start with a fresh LZW dictionary

    for(uint32_t yy=0; yy<m_height; ++yy)
    {
        for(uint32_t xx=0; xx<m_width; ++xx)
        {
            // top-left origin
            uint8_t nextValue = image[(yy*m_width+xx)*4+3];

            // "loser mode" - no compression, every single code is followed immediately by a clear
            //WriteCode( f, stat, nextValue, codeSize );
            //WriteCode( f, stat, 256, codeSize );

            if( curCode < 0 )
            {
                // first value in a new run
                curCode = nextValue;
            }
            else if( codetree[curCode].m_next[nextValue] )
            {
                // current run already in the dictionary
                curCode = codetree[curCode].m_next[nextValue];
            }
            else
            {
                // finish the current run, write a code
                WriteCode(stat, (uint32_t)curCode, codeSize);

                // insert the new run into the dictionary
                codetree[curCode].m_next[nextValue] = (uint16_t)++maxCode;

                if( maxCode >= (1ul << codeSize) )
                {
                    // dictionary entry count has broken a size barrier,
                    // we need more bits for codes
                    codeSize++;
                }
                if( maxCode == 4095 )
                {
                    // the dictionary is full, clear it out and begin anew
                    WriteCode(stat, clearCode, codeSize); // clear tree

                    memset(codetree, 0, sizeof(LzwNode)*4096);
                    codeSize = (uint32_t)(minCodeSize + 1);
                    maxCode = clearCode+1;
                }

                curCode = nextValue;
            }
        }
    }

    // compression footer
    WriteCode(stat, (uint32_t)curCode, codeSize);
    WriteCode(stat, clearCode, codeSize);
    WriteCode(stat, clearCode + 1, (uint32_t)minCodeSize + 1);

    // write out the last partial chunk
    while ( stat.bitIndex ) WriteBit(stat, 0);
    if ( stat.chunkIndex ) WriteChunk(stat);

    fputc(0, m_fpGif); // image block terminator

    free(codetree);
}

int GifGenerator::PickChangedPixels( const uint8_t* lastFrame, uint8_t* frame, int numPixels )
{
    int numChanged = 0;
    uint8_t* writeIter = frame;

    for (int index=0; index<numPixels; ++index)
    {
        if (lastFrame[0] != frame[0] ||
           lastFrame[1] != frame[1] ||
           lastFrame[2] != frame[2])
        {
            writeIter[0] = frame[0];
            writeIter[1] = frame[1];
            writeIter[2] = frame[2];
            ++numChanged;
            writeIter += 4;
        }
        lastFrame += 4;
        frame += 4;
    }

    return numChanged;
}

// Builds a palette by creating a balanced k-d tree of all pixels in the image
void GifGenerator::SplitPalette(uint8_t* image, int numPixels, int firstElt, int lastElt, int splitElt, int splitDist, int treeNode, bool buildForDither, Palette* pal)
{
    if (lastElt <= firstElt || numPixels == 0)
        return;

    // base case, bottom of the tree
    if (lastElt == firstElt+1) {
        if (buildForDither) {
            // Dithering needs at least one color as dark as anything
            // in the image and at least one brightest color -
            // otherwise it builds up error and produces strange artifacts
            if (firstElt == 1) {
                // special case: the darkest color in the image
                int r=255, g=255, b=255;
                for(int index=0; index<numPixels; ++index)
                {
                    r = IMin(r, image[index * 4 + 0]);
                    g = IMin(g, image[index * 4 + 1]);
                    b = IMin(b, image[index * 4 + 2]);
                }

                pal->r[firstElt] = (uint8_t)r;
                pal->g[firstElt] = (uint8_t)g;
                pal->b[firstElt] = (uint8_t)b;

                return;
            }

            if ( firstElt == (1 << pal->bitDepth)-1 ) {
                // special case: the lightest color in the image
                int r=0, g=0, b=0;
                for (int index=0; index<numPixels; ++index) {
                    r = IMax(r, image[index * 4 + 0]);
                    g = IMax(g, image[index * 4 + 1]);
                    b = IMax(b, image[index * 4 + 2]);
                }

                pal->r[firstElt] = (uint8_t)r;
                pal->g[firstElt] = (uint8_t)g;
                pal->b[firstElt] = (uint8_t)b;

                return;
            }
        }

        // otherwise, take the average of all colors in this subcube
        uint64_t r=0, g=0, b=0;
        for (int index=0; index<numPixels; ++index)
        {
            r += image[index*4+0];
            g += image[index*4+1];
            b += image[index*4+2];
        }

        r += (uint64_t)numPixels / 2;  // round to nearest
        g += (uint64_t)numPixels / 2;
        b += (uint64_t)numPixels / 2;

        r /= (uint64_t)numPixels;
        g /= (uint64_t)numPixels;
        b /= (uint64_t)numPixels;

        pal->r[firstElt] = (uint8_t)r;
        pal->g[firstElt] = (uint8_t)g;
        pal->b[firstElt] = (uint8_t)b;

        return;
    }

    // Find the axis with the largest range
    int minR = 255, maxR = 0;
    int minG = 255, maxG = 0;
    int minB = 255, maxB = 0;
    for (int index=0; index<numPixels; ++index)
    {
        int r = image[index*4+0];
        int g = image[index*4+1];
        int b = image[index*4+2];

        if(r > maxR) maxR = r;
        if(r < minR) minR = r;

        if(g > maxG) maxG = g;
        if(g < minG) minG = g;

        if(b > maxB) maxB = b;
        if(b < minB) minB = b;
    }

    int rRange = maxR - minR;
    int gRange = maxG - minG;
    int bRange = maxB - minB;

    // and split along that axis. (incidentally, this means this isn't a "proper" k-d tree but I don't know what else to call it)
    int splitCom = 1;
    if (bRange > gRange) splitCom = 2;
    if (rRange > bRange && rRange > gRange) splitCom = 0;

    int subPixelsA = numPixels * (splitElt - firstElt) / (lastElt - firstElt);
    int subPixelsB = numPixels-subPixelsA;

    PartitionByMedian(image, 0, numPixels, splitCom, subPixelsA);

    pal->treeSplitElt[treeNode] = (uint8_t)splitCom;
    pal->treeSplit[treeNode] = image[subPixelsA*4+splitCom];

    SplitPalette(image, subPixelsA, firstElt, splitElt, splitElt-splitDist, splitDist/2, treeNode*2, buildForDither, pal);
    SplitPalette(image+subPixelsA*4, subPixelsB, splitElt, lastElt,  splitElt+splitDist, splitDist/2, treeNode*2+1, buildForDither, pal);
}

void GifGenerator::GetClosestPaletteColor(Palette* pPal, int r, int g, int b, int& bestInd, int& bestDiff, int treeRoot)
{
    // base case, reached the bottom of the tree
    if (treeRoot > (1<<pPal->bitDepth)-1) {
        int ind = treeRoot-(1<<pPal->bitDepth);
        if (ind == m_kGifTransIndex) return;

        // check whether this color is better than the current winner
        int r_err = r - ((int)pPal->r[ind]);
        int g_err = g - ((int)pPal->g[ind]);
        int b_err = b - ((int)pPal->b[ind]);
        int diff = IAbs(r_err)+IAbs(g_err)+IAbs(b_err);

        if (diff < bestDiff) {
            bestInd = ind;
            bestDiff = diff;
        }

        return;
    }

    // take the appropriate color (r, g, or b) for this node of the k-d tree
    int comps[3]; comps[0] = r; comps[1] = g; comps[2] = b;
    int splitComp = comps[pPal->treeSplitElt[treeRoot]];

    int splitPos = pPal->treeSplit[treeRoot];
    if (splitPos > splitComp) {
        // check the left subtree
        GetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2);
        if ( bestDiff > splitPos - splitComp ) {
            // cannot prove there's not a better value in the right subtree, check that too
            GetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2+1);
        }
    }
    else {
        GetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2+1);
        if ( bestDiff > splitComp - splitPos ) {
            GetClosestPaletteColor(pPal, r, g, b, bestInd, bestDiff, treeRoot*2);
        }
    }
}

void GifGenerator::SwapPixels(uint8_t* image, int pixA, int pixB)
{
    uint8_t rA = image[pixA*4];
    uint8_t gA = image[pixA*4+1];
    uint8_t bA = image[pixA*4+2];
    uint8_t aA = image[pixA*4+3];

    uint8_t rB = image[pixB*4];
    uint8_t gB = image[pixB*4+1];
    uint8_t bB = image[pixB*4+2];
    uint8_t aB = image[pixA*4+3];

    image[pixA*4] = rB;
    image[pixA*4+1] = gB;
    image[pixA*4+2] = bB;
    image[pixA*4+3] = aB;

    image[pixB*4] = rA;
    image[pixB*4+1] = gA;
    image[pixB*4+2] = bA;
    image[pixB*4+3] = aA;
}

// just the partition operation from quicksort
int GifGenerator::Partition(uint8_t* image, const int left, const int right, const int elt, int pivotIndex)
{
    const int pivotValue = image[(pivotIndex)*4+elt];
    SwapPixels(image, pivotIndex, right-1);
    int storeIndex = left;
    bool split = 0;
    for(int index=left; index<right-1; ++index)
    {
        int arrayVal = image[index*4+elt];
        if ( arrayVal < pivotValue ) {
            SwapPixels(image, index, storeIndex);
            ++storeIndex;
        }
        else if( arrayVal == pivotValue ) {
            if (split) {
                SwapPixels(image, index, storeIndex);
                ++storeIndex;
            }
            split = !split;
        }
    }
    SwapPixels(image, storeIndex, right-1);
    return storeIndex;
}

// Perform an incomplete sort, finding all elements above and below the desired median
void GifGenerator::PartitionByMedian(uint8_t* image, int left, int right, int com, int neededCenter)
{
    if (left < right-1) {
        int pivotIndex = left + (right-left)/2;

        pivotIndex = Partition(image, left, right, com, pivotIndex);

        // Only "sort" the section of the array that contains the median
        if (pivotIndex > neededCenter)
            PartitionByMedian(image, left, pivotIndex, com, neededCenter);

        if (pivotIndex < neededCenter)
            PartitionByMedian(image, pivotIndex+1, right, com, neededCenter);
    }
}

// insert a single bit
void GifGenerator::WriteBit( BitStatus& stat, int bit )
{
    bit = bit & 1;
    bit = bit << stat.bitIndex;
    stat.byte |= bit;

    ++stat.bitIndex;
    if ( stat.bitIndex > 7 ) {
        // move the newly-finished byte to the chunk buffer
        stat.chunk[stat.chunkIndex++] = stat.byte;
        // and start a new byte
        stat.bitIndex = 0;
        stat.byte = 0;
    }
}

// write all bytes so far to the file
void GifGenerator::WriteChunk(BitStatus& stat )
{
    fputc((int)stat.chunkIndex, m_fpGif);
    fwrite(stat.chunk, 1, stat.chunkIndex, m_fpGif);

    stat.bitIndex = 0;
    stat.byte = 0;
    stat.chunkIndex = 0;
}

void GifGenerator::WriteCode( BitStatus& stat, int code, int length )
{
    for ( int index=0; index<length; ++index )
    {
        WriteBit(stat, code);
        code = code >> 1;

        if ( stat.chunkIndex == 255 ) {
            WriteChunk(stat);
        }
    }
}



// write a 256-color (8-bit) image palette to the file
void GifGenerator::WritePalette( const Palette* pPal)
{
    fputc(0, m_fpGif);  // first color: transparency
    fputc(0, m_fpGif);
    fputc(0, m_fpGif);

    for(int index=1; index<(1 << pPal->bitDepth); ++index)
    {
        int r = pPal->r[index];
        int g = pPal->g[index];
        int b = pPal->b[index];

        fputc(r, m_fpGif);
        fputc(g, m_fpGif);
        fputc(b, m_fpGif);
    }
}

