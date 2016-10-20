#include <iostream>
#include <fstream>
#include "processors.h"
#include "libfreenect2_data_structures.h"
#include "read_file_handler.h"
#include "device_parameters.h"

int main()
{
    //Load device specific parameters
    DeviceParametersHandler dev;
    dev.init("../parameters/calib_pose_fixed_ir_2.h5");  

    //Load data file containig raw kinect v2 measurents
    ReadFileHandler myfile("../data/bibliotek_view22.h5");
    unsigned char* p0_table;
    int len_p0_table;
    myfile.ReadBuffer(&p0_table, &len_p0_table, "/P0Tables");

    //initialize libreenect2 decoder
    CpuDepthPacketProcessor libfreenect2;
    libfreenect2.loadP0TablesFromCommandResponse(p0_table, len_p0_table);
    libfreenect2.loadXZTables(&dev.xtable[0], &dev.ztable[0]);
    libfreenect2.loadLookupTable(&dev.lut[0]);

    //initialize kde decoder
    CpuKdeDepthPacketProcessor base_kde;
    base_kde.loadP0TablesFromCommandResponse(p0_table, len_p0_table);
    base_kde.loadXZTables(&dev.xtable[0], &dev.ztable[0]);
    base_kde.loadLookupTable(&dev.lut[0]);

    std::fstream out_file_libfreenect2("libfreenect2_depth.bin", std::ios::out | std::ios::binary | std::ios::app);
    std::fstream out_file_libfreenect2_conf("libfreenect2_conf.bin", std::ios::out | std::ios::binary | std::ios::app);
    std::fstream out_file_kde("base_kde_depth.bin", std::ios::out | std::ios::binary | std::ios::app);
    std::fstream out_file_kde_conf("base_kde_conf.bin", std::ios::out | std::ios::binary | std::ios::app);

    //decode and save frames to file
    int frame_num;
    int frame_num_offset = 100;
    int num_frames = 10;
    for(frame_num = frame_num_offset; frame_num < num_frames+frame_num_offset; frame_num++)
    {
        std::cout<<"Frame: "<<frame_num<<std::endl;
        unsigned char *buffer;
        float *depth_base_kde, *depth_libfreenect2, *conf_libfreenect2, *conf_base_kde;
        int len;
        myfile.ReadIrBuffer(&buffer, &len, frame_num);
        
        libfreenect2.process(buffer, &depth_libfreenect2, &conf_libfreenect2);
        base_kde.process(buffer, &depth_base_kde, &conf_base_kde);

        //save buffers to file
        out_file_libfreenect2.write((char*)depth_libfreenect2, 512*424*sizeof(float));
        out_file_libfreenect2_conf.write((char*)conf_libfreenect2, 512*424*sizeof(float));
        out_file_kde.write((char*)depth_base_kde, 512*424*sizeof(float));        
        out_file_kde_conf.write((char*)conf_base_kde, 512*424*sizeof(float));

        delete[] buffer;
        delete[] depth_base_kde;
        delete[] conf_base_kde;   
        delete[] depth_libfreenect2;
        delete[] conf_libfreenect2;
    }

    std::cout<<"Done!\n";
    return 0;
}
