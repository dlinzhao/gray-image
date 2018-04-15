#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <cfloat>
#include <cmath>
#include <cassert>

#include "zlmat.h"

namespace zlmat{

    const int MAX_THICKNESS = 32767;
    typedef unsigned char uchar;

    //! type of line
    enum LineTypes {
        FILLED  = -1,
        LINE_4  = 4, //!< 4-connected line
        LINE_8  = 8, //!< 8-connected line
        LINE_AA = 16
    };

    enum { XY_SHIFT = 16, XY_ONE = 1 << XY_SHIFT, DRAWING_STORAGE_BLOCK = (1<<12) - 256 };

    /* helper macros: filling horizontal row */
    #define ICV_HLINE( ptr, xl, xr, color, pix_size )            \
    {                                                            \
        uchar* hline_ptr = (uchar*)(ptr) + (xl)*(pix_size);      \
        uchar* hline_max_ptr = (uchar*)(ptr) + (xr)*(pix_size);  \
                                                                 \
        for( ; hline_ptr <= hline_max_ptr; hline_ptr += (pix_size))\
        {                                                        \
            int hline_j;                                         \
            for( hline_j = 0; hline_j < (pix_size); hline_j++ )  \
            {                                                    \
                hline_ptr[hline_j] = ((uchar*)color)[hline_j];   \
            }                                                    \
        }                                                        \
    }

    bool clipLine(Size img_size, Point& pt1, Point& pt2);
    bool clipLine(Size2l img_size, Point2l& pt1, Point2l& pt2);
    bool clipLine(Rect img_rect, Point& pt1, Point& pt2);

    void ThickLine( Mat& img, Point2l p0, Point2l p1, const void* color,
               int thickness, int line_type, int flags, int shift );
    /* filling convex polygon. v - array of vertices, ntps - number of points */
    void FillConvexPoly( Mat& img, const Point2l* v, int npts,
                         const void* color, int line_type, int shift );
    /* draws simple or filled circle */
    void Circle( Mat& img, Point center, int radius, const void* color, int fill );
    void Line( Mat& img, Point pt1, Point pt2, const void* _color, int connectivity = 8 );
    void Line2( Mat& img, Point2l pt1, Point2l pt2, const void* color);

    bool clipLine(Size img_size, Point& pt1, Point& pt2)
    {
        Point2l p1(pt1);
        Point2l p2(pt2);
        bool inside = clipLine(Size2l((long long)img_size.width,
                                         (long long)img_size.height), p1, p2);
        pt1.x = (int)p1.x;
        pt1.y = (int)p1.y;
        pt2.x = (int)p2.x;
        pt2.y = (int)p2.y;
        return inside;
    }

    bool clipLine( Size2l img_size, Point2l& pt1, Point2l& pt2)
    {
        int c1, c2;
        long long right = img_size.width-1, bottom = img_size.height-1;

        if( img_size.width <= 0 || img_size.height <= 0 )
            return false;

        long long &x1 = pt1.x, &y1 = pt1.y, &x2 = pt2.x, &y2 = pt2.y;
        c1 = (x1 < 0) + (x1 > right) * 2 + (y1 < 0) * 4 + (y1 > bottom) * 8;
        c2 = (x2 < 0) + (x2 > right) * 2 + (y2 < 0) * 4 + (y2 > bottom) * 8;

        if( (c1 & c2) == 0 && (c1 | c2) != 0 )
        {
            long long a;
            if( c1 & 12 )
            {
                a = c1 < 8 ? 0 : bottom;
                x1 +=  (a - y1) * (x2 - x1) / (y2 - y1);
                y1 = a;
                c1 = (x1 < 0) + (x1 > right) * 2;
            }
            if( c2 & 12 )
            {
                a = c2 < 8 ? 0 : bottom;
                x2 += (a - y2) * (x2 - x1) / (y2 - y1);
                y2 = a;
                c2 = (x2 < 0) + (x2 > right) * 2;
            }
            if( (c1 & c2) == 0 && (c1 | c2) != 0 )
            {
                if( c1 )
                {
                    a = c1 == 1 ? 0 : right;
                    y1 += (a - x1) * (y2 - y1) / (x2 - x1);
                    x1 = a;
                    c1 = 0;
                }
                if( c2 )
                {
                    a = c2 == 1 ? 0 : right;
                    y2 += (a - x2) * (y2 - y1) / (x2 - x1);
                    x2 = a;
                    c2 = 0;
                }
            }

            assert( (c1 & c2) != 0 || (x1 | y1 | x2 | y2) >= 0 );
        }

        return (c1 | c2) == 0;
    }

    bool clipLine(Rect img_rect, Point& pt1, Point& pt2 )
    {
        Point tl = img_rect.tl();
        pt1 -= tl; pt2 -= tl;
        bool inside = clipLine(img_rect.size(), pt1, pt2);
        pt1 += tl; pt2 += tl;

        return inside;
    }

    // === Mat implementation ===

    Mat_::Mat_()
    {
        cols = 0;
        rows = 0;
        step = 0;
        data = nullptr;
    }

    Mat_::Mat_(int rows, int cols)
    {
        create(rows, cols);
    }

    Mat_::Mat_(Size size)
    {
        create(size.height, size.width);
    }

    Mat_::Mat_(const Mat_ &m)
    {
        create(m.rows, m.cols);
        memcpy(data, m.data, rows * cols);
    }

    void Mat_::fillValue(unsigned char value)
    {
        std::fill(data, data + rows * cols, value);
    }

    void Mat_::create(int rows, int cols)
    {
        this->rows = rows;
        this->cols = cols;
        this->step = cols;
        this->data = new unsigned char[rows * cols];
    }

    void Mat_::resize(Mat_ &dst, Size s)
    {
        double scale_x = (double)cols / s.width;
        double scale_y = (double)rows / s.height;

        if(dst.isEmpty())
            dst.create(s.height, s.width);

        unsigned char* pdst = dst.data;
        int dstep = dst.step;

        for (int j = 0; j < dst.rows; ++j)
        {
            float fy = (float)((j + 0.5) * scale_y - 0.5);
            int sy = floor(fy);
            fy -= sy;
            sy = std::min(sy, rows - 2);
            sy = std::max(0, sy);

            short cbufy[2];
            cbufy[0] = std::max(std::min(int(1.f - fy) * 2048, 2048), 0); //saturate_cast<short>
            cbufy[1] = 2048 - cbufy[0];

            for (int i = 0; i < dst.cols; ++i)
            {
                float fx = (float)((i + 0.5) * scale_x - 0.5);
                int sx = floor(fx);
                fx -= sx;

                if (sx < 0) {
                    fx = 0, sx = 0;
                }
                if (sx >= cols - 1) {
                    fx = 0, sx = cols - 2;
                }

                short cbufx[2];
                cbufx[0] = std::max(std::min(int(1.f - fx) * 2048, 2048), 0); //saturate_cast<short>
                cbufx[1] = 2048 - cbufx[0];

                *(pdst+ j*dstep + i) = (*(data + sy*step + sx) * cbufx[0] * cbufy[0] +
                    *(data + (sy+1)*step + sx) * cbufx[0] * cbufy[1] +
                    *(data + sy*step + (sx+1)) * cbufx[1] * cbufy[0] +
                    *(data + (sy+1)*step + (sx+1)) * cbufx[1] * cbufy[1]) >> 22;
            }
        }
    }

    bool Mat_::isEmpty()
    {
        if(data == nullptr)
            return true;
        return false;
    }

    Mat_::~Mat_()
    {
        if(data != nullptr)
        {
            delete[] data;
            data = nullptr;
        }
    }

    std::ostream & operator << (std::ostream &out, const Mat_ &m)
    {
        out << "Mat [" << m.rows << " x " << m.rows << "]:\n";
        for(int i = 0; i < m.rows; i++)
        {
            for(int j = 0; j < m.cols; j++)
            {
                out << (int)m.data[i * m.step + j] << "\t";
            }
            out << "\n";
        }
        return out;
    }

    /*
       Initializes line iterator.
       Returns number of points on the line or negative number if error.
    */
    LineIterator::LineIterator(Mat& img, Point pt1, Point pt2,
                               int connectivity, bool left_to_right)
    {
        count = -1;

        assert( connectivity == 8 || connectivity == 4 );

        if( (unsigned)pt1.x >= (unsigned)(img.cols) ||
            (unsigned)pt2.x >= (unsigned)(img.cols) ||
            (unsigned)pt1.y >= (unsigned)(img.rows) ||
            (unsigned)pt2.y >= (unsigned)(img.rows) )
        {
            if( !clipLine(img.size(), pt1, pt2 ) )
            {
                ptr = img.data;
                err = plusDelta = minusDelta = plusStep = minusStep = count = 0;
                return;
            }
        }

        int bt_pix0 = (int)img.elemSize(), bt_pix = bt_pix0;
        size_t istep = img.step;

        int dx = pt2.x - pt1.x;
        int dy = pt2.y - pt1.y;
        int s = dx < 0 ? -1 : 0;

        if( left_to_right )
        {
            dx = (dx ^ s) - s;
            dy = (dy ^ s) - s;
            pt1.x ^= (pt1.x ^ pt2.x) & s;
            pt1.y ^= (pt1.y ^ pt2.y) & s;
        }
        else
        {
            dx = (dx ^ s) - s;
            bt_pix = (bt_pix ^ s) - s;
        }

        ptr = (uchar*)(img.data + pt1.y * istep + pt1.x * bt_pix0);

        s = dy < 0 ? -1 : 0;
        dy = (dy ^ s) - s;
        istep = (istep ^ s) - s;

        s = dy > dx ? -1 : 0;

        /* conditional swaps */
        dx ^= dy & s;
        dy ^= dx & s;
        dx ^= dy & s;

        bt_pix ^= istep & s;
        istep ^= bt_pix & s;
        bt_pix ^= istep & s;

        if( connectivity == 8 )
        {
            assert( dx >= 0 && dy >= 0 );

            err = dx - (dy + dy);
            plusDelta = dx + dx;
            minusDelta = -(dy + dy);
            plusStep = (int)istep;
            minusStep = bt_pix;
            count = dx + 1;
        }
        else /* connectivity == 4 */
        {
            assert( dx >= 0 && dy >= 0 );

            err = 0;
            plusDelta = (dx + dx) + (dy + dy);
            minusDelta = -(dy + dy);
            plusStep = (int)istep - bt_pix;
            minusStep = bt_pix;
            count = dx + dy + 1;
        }

        this->ptr0 = img.ptr();
        this->step = (int)img.step;
        this->elemSize = bt_pix0;
    }


    void line(Mat &_img, Point pt1, Point pt2, const Color color,
               int thickness, int line_type, int shift )
    {
        assert( 0 <= thickness && thickness <= MAX_THICKNESS );
        assert( 0 <= shift && shift <= XY_SHIFT );

        unsigned char buf[4] = {color.c};

        ThickLine( _img, pt1, pt2, buf, thickness, line_type, 3, shift );
    }

    void ThickLine( Mat& img, Point2l p0, Point2l p1, const void* color,
               int thickness, int line_type, int flags, int shift )
    {
        static const double INV_XY_ONE = 1./XY_ONE;

        p0.x <<= XY_SHIFT - shift;
        p0.y <<= XY_SHIFT - shift;
        p1.x <<= XY_SHIFT - shift;
        p1.y <<= XY_SHIFT - shift;

        Point2l pt[4], dp = Point2l(0,0);
        double dx = (p0.x - p1.x)*INV_XY_ONE, dy = (p1.y - p0.y)*INV_XY_ONE;
        double r = dx * dx + dy * dy;
        int i, oddThickness = thickness & 1;
        thickness <<= XY_SHIFT - 1;

        if( fabs(r) > DBL_EPSILON )
        {
            r = (thickness + oddThickness*XY_ONE*0.5)/std::sqrt(r);
            dp.x = round( dy * r );
            dp.y = round( dx * r );

            pt[0].x = p0.x + dp.x;
            pt[0].y = p0.y + dp.y;
            pt[1].x = p0.x - dp.x;
            pt[1].y = p0.y - dp.y;
            pt[2].x = p1.x - dp.x;
            pt[2].y = p1.y - dp.y;
            pt[3].x = p1.x + dp.x;
            pt[3].y = p1.y + dp.y;

            FillConvexPoly( img, pt, 4, color, line_type, XY_SHIFT );
        }

        for( i = 0; i < 2; i++ )
        {
            if( flags & (i+1) )
            {
                Point center;
                center.x = (int)((p0.x + (XY_ONE>>1)) >> XY_SHIFT);
                center.y = (int)((p0.y + (XY_ONE>>1)) >> XY_SHIFT);
                Circle( img, center, (thickness + (XY_ONE>>1)) >> XY_SHIFT, color, 1 );
            }
            p0 = p1;
        }
    }

    /* filling convex polygon. v - array of vertices, ntps - number of points */
    void FillConvexPoly( Mat& img, const Point2l* v,
                         int npts, const void* color, int line_type, int shift )
    {
        struct
        {
            int idx, di;
            long long x, dx;
            int ye;
        }edge[2];

        int delta = 1 << shift >> 1;
        int i, y, imin = 0;
        int edges = npts;
        long long xmin, xmax, ymin, ymax;
        uchar* ptr = img.ptr();
        Size size = img.size();
        int pix_size = (int)img.elemSize();
        Point2l p0;
        int delta1, delta2;

        delta1 = delta2 = XY_ONE >> 1;

        p0 = v[npts - 1];
        p0.x <<= XY_SHIFT - shift;
        p0.y <<= XY_SHIFT - shift;

        assert( 0 <= shift && shift <= XY_SHIFT );
        xmin = xmax = v[0].x;
        ymin = ymax = v[0].y;

        for( i = 0; i < npts; i++ )
        {
            Point2l p = v[i];
            if( p.y < ymin )
            {
                ymin = p.y;
                imin = i;
            }

            ymax = std::max( ymax, p.y );
            xmax = std::max( xmax, p.x );
            xmin = MIN( xmin, p.x );

            p.x <<= XY_SHIFT - shift;
            p.y <<= XY_SHIFT - shift;

            if( line_type <= 8 )
            {
                if( shift == 0 )
                {
                    Point pt0, pt1;
                    pt0.x = (int)(p0.x >> XY_SHIFT);
                    pt0.y = (int)(p0.y >> XY_SHIFT);
                    pt1.x = (int)(p.x >> XY_SHIFT);
                    pt1.y = (int)(p.y >> XY_SHIFT);
                    Line( img, pt0, pt1, color, line_type );
                }
                else
                    Line2( img, p0, p, color );
            }
            p0 = p;
        }

        xmin = (xmin + delta) >> shift;
        xmax = (xmax + delta) >> shift;
        ymin = (ymin + delta) >> shift;
        ymax = (ymax + delta) >> shift;

        if( npts < 3 || (int)xmax < 0 || (int)ymax < 0 || (int)xmin >= size.width || (int)ymin >= size.height )
            return;

        ymax = MIN( ymax, size.height - 1 );
        edge[0].idx = edge[1].idx = imin;

        edge[0].ye = edge[1].ye = y = (int)ymin;
        edge[0].di = 1;
        edge[1].di = npts - 1;

        ptr += img.step*y;

        do
        {
            if( line_type < LINE_AA || y < (int)ymax || y == (int)ymin )
            {
                for( i = 0; i < 2; i++ )
                {
                    if( y >= edge[i].ye )
                    {
                        int idx0 = edge[i].idx, di = edge[i].di;
                        int idx = idx0 + di;
                        if (idx >= npts) idx -= npts;
                        int ty = 0;

                        for (; edges-- > 0; )
                        {
                            ty = (int)((v[idx].y + delta) >> shift);
                            if (ty > y)
                            {
                                long long xs = v[idx0].x;
                                long long xe = v[idx].x;
                                if (shift != XY_SHIFT)
                                {
                                    xs <<= XY_SHIFT - shift;
                                    xe <<= XY_SHIFT - shift;
                                }

                                edge[i].ye = ty;
                                edge[i].dx = ((xe - xs)*2 + (ty - y)) / (2 * (ty - y));
                                edge[i].x = xs;
                                edge[i].idx = idx;
                                break;
                            }
                            idx0 = idx;
                            idx += di;
                            if (idx >= npts) idx -= npts;
                        }
                    }
                }
            }

            if (y >= 0)
            {
                int left = 0, right = 1;
                if (edge[0].x > edge[1].x)
                {
                    left = 1, right = 0;
                }

                int xx1 = (int)((edge[left].x + delta1) >> XY_SHIFT);
                int xx2 = (int)((edge[right].x + delta2) >> XY_SHIFT);

                if( xx2 >= 0 && xx1 < size.width )
                {
                    if( xx1 < 0 )
                        xx1 = 0;
                    if( xx2 >= size.width )
                        xx2 = size.width - 1;
                    ICV_HLINE( ptr, xx1, xx2, color, pix_size );
                }
            }
            else
            {
                // TODO optimize scan for negative y
            }

            edge[0].x += edge[0].dx;
            edge[1].x += edge[1].dx;
            ptr += img.step;
        }
        while( ++y <= (int)ymax );
    }
    /* draws simple or filled circle */
    void Circle( Mat& img, Point center, int radius, const void* color, int fill )
    {
        Size size = img.size();
        size_t step = img.step;
        int pix_size = (int)img.elemSize();
        uchar* ptr = img.ptr();
        int err = 0, dx = radius, dy = 0, plus = 1, minus = (radius << 1) - 1;
        int inside = center.x >= radius && center.x < size.width - radius &&
            center.y >= radius && center.y < size.height - radius;

        #define ICV_PUT_POINT( ptr, x )     \
            memcpy( ptr + (x)*pix_size, color, pix_size );

        while( dx >= dy )
        {
            int mask;
            int y11 = center.y - dy, y12 = center.y + dy, y21 = center.y - dx, y22 = center.y + dx;
            int x11 = center.x - dx, x12 = center.x + dx, x21 = center.x - dy, x22 = center.x + dy;

            if( inside )
            {
                uchar *tptr0 = ptr + y11 * step;
                uchar *tptr1 = ptr + y12 * step;

                if( !fill )
                {
                    ICV_PUT_POINT( tptr0, x11 );
                    ICV_PUT_POINT( tptr1, x11 );
                    ICV_PUT_POINT( tptr0, x12 );
                    ICV_PUT_POINT( tptr1, x12 );
                }
                else
                {
                    ICV_HLINE( tptr0, x11, x12, color, pix_size );
                    ICV_HLINE( tptr1, x11, x12, color, pix_size );
                }

                tptr0 = ptr + y21 * step;
                tptr1 = ptr + y22 * step;

                if( !fill )
                {
                    ICV_PUT_POINT( tptr0, x21 );
                    ICV_PUT_POINT( tptr1, x21 );
                    ICV_PUT_POINT( tptr0, x22 );
                    ICV_PUT_POINT( tptr1, x22 );
                }
                else
                {
                    ICV_HLINE( tptr0, x21, x22, color, pix_size );
                    ICV_HLINE( tptr1, x21, x22, color, pix_size );
                }
            }
            else if( x11 < size.width && x12 >= 0 && y21 < size.height && y22 >= 0 )
            {
                if( fill )
                {
                    x11 = std::max( x11, 0 );
                    x12 = MIN( x12, size.width - 1 );
                }

                if( (unsigned)y11 < (unsigned)size.height )
                {
                    uchar *tptr = ptr + y11 * step;

                    if( !fill )
                    {
                        if( x11 >= 0 )
                            ICV_PUT_POINT( tptr, x11 );
                        if( x12 < size.width )
                            ICV_PUT_POINT( tptr, x12 );
                    }
                    else
                        ICV_HLINE( tptr, x11, x12, color, pix_size );
                }

                if( (unsigned)y12 < (unsigned)size.height )
                {
                    uchar *tptr = ptr + y12 * step;

                    if( !fill )
                    {
                        if( x11 >= 0 )
                            ICV_PUT_POINT( tptr, x11 );
                        if( x12 < size.width )
                            ICV_PUT_POINT( tptr, x12 );
                    }
                    else
                        ICV_HLINE( tptr, x11, x12, color, pix_size );
                }

                if( x21 < size.width && x22 >= 0 )
                {
                    if( fill )
                    {
                        x21 = std::max( x21, 0 );
                        x22 = MIN( x22, size.width - 1 );
                    }

                    if( (unsigned)y21 < (unsigned)size.height )
                    {
                        uchar *tptr = ptr + y21 * step;

                        if( !fill )
                        {
                            if( x21 >= 0 )
                                ICV_PUT_POINT( tptr, x21 );
                            if( x22 < size.width )
                                ICV_PUT_POINT( tptr, x22 );
                        }
                        else
                            ICV_HLINE( tptr, x21, x22, color, pix_size );
                    }

                    if( (unsigned)y22 < (unsigned)size.height )
                    {
                        uchar *tptr = ptr + y22 * step;

                        if( !fill )
                        {
                            if( x21 >= 0 )
                                ICV_PUT_POINT( tptr, x21 );
                            if( x22 < size.width )
                                ICV_PUT_POINT( tptr, x22 );
                        }
                        else
                            ICV_HLINE( tptr, x21, x22, color, pix_size );
                    }
                }
            }
            dy++;
            err += plus;
            plus += 2;

            mask = (err <= 0) - 1;

            err -= minus & mask;
            dx += mask;
            minus -= mask & 2;
        }

        #undef  ICV_PUT_POINT
    }

    void Line( Mat& img, Point pt1, Point pt2, const void* _color, int connectivity)
    {
        LineIterator iterator(img, pt1, pt2, connectivity, true);
        int i, count = iterator.count;
        int pix_size = (int)img.elemSize();
        const uchar* color = (const uchar*)_color;

        for( i = 0; i < count; i++, ++iterator )
        {
            uchar* ptr = *iterator;
            if( pix_size == 1 )
                ptr[0] = color[0];
            else if( pix_size == 3 )
            {
                ptr[0] = color[0];
                ptr[1] = color[1];
                ptr[2] = color[2];
            }
            else
                memcpy( *iterator, color, pix_size );
        }
    }

    void  Line2( Mat& img, Point2l pt1, Point2l pt2, const void* color)
    {
        long long dx, dy;
        int ecount;
        long long ax, ay;
        long long i, j;
        int x, y;
        long long x_step, y_step;
        int cb = ((uchar*)color)[0];
        int cg = ((uchar*)color)[1];
        int cr = ((uchar*)color)[2];
        int pix_size = (int)img.elemSize();
        uchar *ptr = img.ptr(), *tptr;
        size_t step = img.step;
        Size size = img.size();

        long long sw = ((long long)size.width) << XY_SHIFT;
        long long sh = ((long long)size.height) << XY_SHIFT;
        Size2l sizeScaled(sw, sh);

        if(clipLine(sizeScaled, pt1, pt2) == false)
            return;

        dx = pt2.x - pt1.x;
        dy = pt2.y - pt1.y;

        j = dx < 0 ? -1 : 0;
        ax = (dx ^ j) - j;
        i = dy < 0 ? -1 : 0;
        ay = (dy ^ i) - i;

        if( ax > ay )
        {
            dx = ax;
            dy = (dy ^ j) - j;
            pt1.x ^= pt2.x & j;
            pt2.x ^= pt1.x & j;
            pt1.x ^= pt2.x & j;
            pt1.y ^= pt2.y & j;
            pt2.y ^= pt1.y & j;
            pt1.y ^= pt2.y & j;

            x_step = XY_ONE;
            y_step = (dy << XY_SHIFT) / (ax | 1);
            ecount = (int)((pt2.x - pt1.x) >> XY_SHIFT);
        }
        else
        {
            dy = ay;
            dx = (dx ^ i) - i;
            pt1.x ^= pt2.x & i;
            pt2.x ^= pt1.x & i;
            pt1.x ^= pt2.x & i;
            pt1.y ^= pt2.y & i;
            pt2.y ^= pt1.y & i;
            pt1.y ^= pt2.y & i;

            x_step = (dx << XY_SHIFT) / (ay | 1);
            y_step = XY_ONE;
            ecount = (int)((pt2.y - pt1.y) >> XY_SHIFT);
        }

        pt1.x += (XY_ONE >> 1);
        pt1.y += (XY_ONE >> 1);

        if( pix_size == 3 )
        {
            #define  ICV_PUT_POINT(_x,_y)   \
            x = (_x); y = (_y);             \
            if( 0 <= x && x < size.width && \
                0 <= y && y < size.height ) \
            {                               \
                tptr = ptr + y*step + x*3;  \
                tptr[0] = (uchar)cb;        \
                tptr[1] = (uchar)cg;        \
                tptr[2] = (uchar)cr;        \
            }

            ICV_PUT_POINT((int)((pt2.x + (XY_ONE >> 1)) >> XY_SHIFT),
                          (int)((pt2.y + (XY_ONE >> 1)) >> XY_SHIFT));

            if( ax > ay )
            {
                pt1.x >>= XY_SHIFT;

                while( ecount >= 0 )
                {
                    ICV_PUT_POINT((int)(pt1.x), (int)(pt1.y >> XY_SHIFT));
                    pt1.x++;
                    pt1.y += y_step;
                    ecount--;
                }
            }
            else
            {
                pt1.y >>= XY_SHIFT;

                while( ecount >= 0 )
                {
                    ICV_PUT_POINT((int)(pt1.x >> XY_SHIFT), (int)(pt1.y));
                    pt1.x += x_step;
                    pt1.y++;
                    ecount--;
                }
            }

            #undef ICV_PUT_POINT
        }
        else if( pix_size == 1 )
        {
            #define  ICV_PUT_POINT(_x,_y) \
            x = (_x); y = (_y);           \
            if( 0 <= x && x < size.width && \
                0 <= y && y < size.height ) \
            {                           \
                tptr = ptr + y*step + x;\
                tptr[0] = (uchar)cb;    \
            }

            ICV_PUT_POINT((int)((pt2.x + (XY_ONE >> 1)) >> XY_SHIFT),
                          (int)((pt2.y + (XY_ONE >> 1)) >> XY_SHIFT));

            if( ax > ay )
            {
                pt1.x >>= XY_SHIFT;

                while( ecount >= 0 )
                {
                    ICV_PUT_POINT((int)(pt1.x), (int)(pt1.y >> XY_SHIFT));
                    pt1.x++;
                    pt1.y += y_step;
                    ecount--;
                }
            }
            else
            {
                pt1.y >>= XY_SHIFT;

                while( ecount >= 0 )
                {
                    ICV_PUT_POINT((int)(pt1.x >> XY_SHIFT), (int)(pt1.y));
                    pt1.x += x_step;
                    pt1.y++;
                    ecount--;
                }
            }

            #undef ICV_PUT_POINT
        }
        else
        {
            #define  ICV_PUT_POINT(_x,_y)   \
            x = (_x); y = (_y);             \
            if( 0 <= x && x < size.width && \
                0 <= y && y < size.height ) \
            {                               \
                tptr = ptr + y*step + x*pix_size;\
                for( j = 0; j < pix_size; j++ ) \
                    tptr[j] = ((uchar*)color)[j]; \
            }

            ICV_PUT_POINT((int)((pt2.x + (XY_ONE >> 1)) >> XY_SHIFT),
                          (int)((pt2.y + (XY_ONE >> 1)) >> XY_SHIFT));

            if( ax > ay )
            {
                pt1.x >>= XY_SHIFT;

                while( ecount >= 0 )
                {
                    ICV_PUT_POINT((int)(pt1.x), (int)(pt1.y >> XY_SHIFT));
                    pt1.x++;
                    pt1.y += y_step;
                    ecount--;
                }
            }
            else
            {
                pt1.y >>= XY_SHIFT;

                while( ecount >= 0 )
                {
                    ICV_PUT_POINT((int)(pt1.x >> XY_SHIFT), (int)(pt1.y));
                    pt1.x += x_step;
                    pt1.y++;
                    ecount--;
                }
            }
            #undef ICV_PUT_POINT
        }
    }
}
