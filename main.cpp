#include "zlmat.h"

#include <ctime>
#include <cstdlib>
#include <iostream>

int main(int argc, char *argv[])
{
    if(argc != 3)
    {
        std::cerr << "Useage: " << argv[0] << " rows cols" << std::endl;
        return -1;
    }
    int rows = atoi(argv[1]);
    int cols = atoi(argv[2]);

    srand(time(0));
    // create a image
    zlmat::Mat img(rows, cols);
    std::cout << img << std::endl;
    // generate some points
    zlmat::Point pt1, pt2;
    pt1.x = rand() % cols;
    pt1.y = rand() % rows;

    for(int i = 0; i < 2; i++)
    {
        pt2.x = rand() % cols;
        pt2.y = rand() % rows;
        zlmat::line(img, pt1, pt2, zlmat::Color(255), 3);
        pt1 = pt2;
    }
    // print img
    std::cout << img << std::endl;
    return 0;
}
