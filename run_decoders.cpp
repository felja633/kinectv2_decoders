#include <iostream>
#include "processors.h"
#include "libfreenect2_data_structures.h"
#include "read_file_handler.h"
#include "device_parameters.h"

int main()
{
    DeviceParametersHandler dev;
    dev.init("../parameters/calib_pose_fixed_ir_2.h5");
    
    std::cout<<"funkar\n";   

    ReadFileHandler myfile("../data/bibliotek_view22.h5");
    unsigned char* p0_table;
    int len_p0_table;
    myfile.ReadBuffer(&p0_table, &len_p0_table, "/P0Tables");

    CpuDepthPacketProcessor libfreenect2;
    libfreenect2.loadP0TablesFromCommandResponse(p0_table, len_p0_table);
    libfreenect2.loadXZTables(dev.ztable, dev.xtable);
    libfreenect2.loadLookupTable(dev.lut);

    CpuKdeDepthPacketProcessor base_kde;
    base_kde.loadP0TablesFromCommandResponse(p0_table, len_p0_table);
    base_kde.loadXZTables(dev.ztable, dev.xtable);
    base_kde.loadLookupTable(dev.lut);

    int frame_num = 0;
    int num_frames = 10;
    for(frame_num = 0; frame_num < num_frames; frame_num++)
    {
        std::cout<<"Frane: "<<frame_num<<std::endl;
        unsigned char *buffer;
        float *depth_base, *depth_libfreenect2, *ir_libfreenect2, *ir_base_kde;
        int len;
        myfile.ReadIrBuffer(&buffer, &len, frame_num);
        
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
