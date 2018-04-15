#ifndef __ZLMAT_H_
#define __ZLMAT_H_

#include <algorithm>
#include <cstring>
#include <ostream>
#include <cstddef>

namespace zlmat {

    #ifndef MIN
    template<typename T1, typename T2>
    inline auto MIN(T1 x, T2 a)-> decltype(x <= a? x: a)
    { return x <= a? x: a; }
    #endif

    template <typename T>
    struct Point_
    {
        T x, y;
        Point_():x(0), y(0){}
        Point_(const Point_ &p): x(p.x), y(p.y){}
        Point_(T _x, T _y):x(_x), y(_y){}

        template<typename T2> operator Point_<T2>() const;
        void operator -= (Point_ p){x -= p.x; y -= p.y;}
        void operator += (Point_ p){x += p.x; y += p.y;}
    };

    typedef Point_<int> Point;
    typedef Point_<long long> Point2l;

    template <typename T>
    class Size_
    {
    public:
        T width, height;
        Size_():width(0), height(0){}
        Size_(T w, T h):width(w), height(h){}

        template<typename T2> operator Size_<T2>() const;
    };

    typedef Size_<int> Size;
    typedef Size_<long long> Size2l;

    template <typename T>
    class Rect_
    {
    public:
        T x, y, width, height;
        Rect_(){}
        Rect_(T _x, T _y, T w, T h):
            x(_x), y(_y), width(w), height(h){}
        //! the top-left corner
        Point_<T> tl() const {return Point_<T>(x, y);}
        Size_<T> size() const {return Size_<T>(width, height);}
    };

    typedef Rect_<int> Rect;

    class Mat_
    {
    public:
        Mat_();
        Mat_(int rows, int cols);
        Mat_(Size size);
        Mat_(const Mat_ &m);
        ~Mat_();

        friend std::ostream & operator << (std::ostream &out, const Mat_ &m);

        void fillValue(unsigned char value=0);
        void create(int rows, int cols);
        void resize(Mat_ &dst, Size s);
        bool isEmpty();

        inline unsigned char * ptr(){return data;}
        Size size(){return Size(cols, rows);}
        size_t elemSize(){return 1;}

    public:
        int cols;
        int rows;
        int step;
        unsigned char *data;
    };

    typedef Mat_ Mat;

    typedef struct Color_
    {
        unsigned char c;
        Color_(){}
        Color_(unsigned char a){c = a;}
    }Color;

    class LineIterator
    {
    public:
        LineIterator( Mat& img, Point pt1, Point pt2,
                      int connectivity = 8, bool leftToRight = false );
        /** @brief returns pointer to the current pixel
        */
        unsigned char *operator *();
        /** @brief prefix increment operator (++it). shifts iterator to the next pixel
        */
        LineIterator& operator ++();
        /** @brief postfix increment operator (it++). shifts iterator to the next pixel
        */
        LineIterator operator ++(int);
        /** @brief returns coordinates of the current pixel
        */
        Point pos() const;

        unsigned char* ptr;
        const unsigned char* ptr0;
        int step, elemSize;
        int err, count;
        int minusDelta, plusDelta;
        int minusStep, plusStep;
    };

    // === Point implementation ===
    template<typename T> template<typename T2> inline
    Point_<T>::operator Point_<T2>() const
    {
        return Point_<T2>(static_cast<T2>(x), static_cast<T2>(y));
    }

    template<typename T> template<typename T2> inline
    Size_<T>::operator Size_<T2>() const
    {
        return Size_<T2>(static_cast<T2>(width), static_cast<T2>(height));
    }

    // === LineIterator implementation ===

    inline
    unsigned char* LineIterator::operator *()
    {
        return ptr;
    }

    inline
    LineIterator& LineIterator::operator ++()
    {
        int mask = err < 0 ? -1 : 0;
        err += minusDelta + (plusDelta & mask);
        ptr += minusStep + (plusStep & mask);
        return *this;
    }

    inline
    LineIterator LineIterator::operator ++(int)
    {
        LineIterator it = *this;
        ++(*this);
        return it;
    }

    inline
    Point LineIterator::pos() const
    {
        Point p;
        p.y = (int)((ptr - ptr0)/step);
        p.x = (int)(((ptr - ptr0) - p.y*step)/elemSize);
        return p;
    }

    void line(Mat &_img, Point pt1, Point pt2, const Color color,
               int thickness=1, int line_type=8, int shift=0 );
}
#endif // LINE_H
