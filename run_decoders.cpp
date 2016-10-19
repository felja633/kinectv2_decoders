#include <iostream>
#include "processors.h"
#include "libfreenect2_data_structures.h"
#include "read_file_handler.h"

int main()
{
    ReadFileHandler myfile("name.h5");
    std::cout<<"hello\n";
    return 0;
}
