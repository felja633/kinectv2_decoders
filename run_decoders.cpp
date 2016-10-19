#include <iostream>
#include "processors.h"
#include "libfreenect2_data_structures.h"
#include "read_file_handler.h"
#include "device_parameters.h"

int main()
{
    ReadFileHandler myfile("name.h5");
    unsigned char* p0_table;
    int len_p0_table;
    myfile.ReadBuffer(&p0_table, &len_p0_table, "p0_table");
    DeviceParametersHandler dev;
    dev.init("parameters/calib_fixed_ir_2.5");

    CpuDepthPacketProcessor libfreenect2();
    libfreenect2.loadP0TablesFromCommandResponse(p0_table, len_p0_table);
    libfreenect2.loadXZTables(dev.ztables, dev.xtables);
    libfreenect2.loadLookupTable(dev.lut);

    CpuKdeDepthPacketProcessor base_kde();
    base_kde.loadP0TablesFromCommandResponse(p0_table, len_p0_table);
    base_kde.loadXZTables(dev.ztables, dev.xtables);
    base_kde.loadLookupTable(dev.lut);

    int frame_num = 0;
    int num_frames = 10;
    for(frame_num = 0; frame_num < num_frames; frame_num++)
    {
        float* buffer, depth_base, depth_libfreenect2, ir_libfreenect2, ir_base_kde;
        int len;
        ReadIrBuffer(&buffer, &len, frame_num);
        
        libfreenect2.process(buffer, &depth_libfreenect2, &ir_libfreenect2);
        base_kde.process(buffer, &depth_base, &ir_base_kde);

        //save buffers to file

        delete[] buffer;
        delete[] depth_base;
        delete[] ir_base_kde;   
        delete[] depth_libfreenect2;
        delete[] ir_libfreenect2;
    }

    std::cout<<"hello\n";
    return 0;
}
