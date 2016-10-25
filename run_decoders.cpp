#include <iostream>
#include <fstream>
#include "processors.h"
#include "libfreenect2_data_structures.h"
#include "read_file_handler.h"
#include "device_parameters.h"
#include "setup_parser.h"
#include <stdio.h>
#include <string.h>

int main(int argc, char** argv)
{
    if(argc != 4)
    {
        std::cout<<"usage:\n ./kinectv2_decoders /path/to/convig/setup.xml dataset /path/to/dataset/";
        return 0;
    }

    std::string setup_xmlfile = argv[1];
    std::string data_set = argv[2];
    std::string data_set_path = argv[3];

    std::string data_file_name = data_set_path;
    data_file_name.append(data_set);
    data_file_name.append(".h5");

    //parse parameter file
    SetupParser parser;
    parser.init(setup_xmlfile);
    
    std::vector<Parameters> params = parser.getParameters();
  
    //Load device specific parameters
    DeviceParametersHandler dev;
    dev.init("../parameters/calib_pose_fixed_ir_2.h5");  

    //Load data file containig raw kinect v2 measurents
    ReadFileHandler myfile(data_file_name);
    unsigned char* p0_table;
    int len_p0_table;
    myfile.ReadBuffer(&p0_table, &len_p0_table, "/P0Tables");

    std::vector<DepthPacketProcessor*> processors;
    std::vector<std::fstream*> depth_files, conf_files; 
    for(unsigned int i = 0; i < params.size(); ++i)
    {
        if(strcmp(params[i].pipeline.c_str(),"libfreenect2") == 0)
        {
            //initialize libreenect2 decoder
            CpuDepthPacketProcessor* libfreenect2 = new CpuDepthPacketProcessor();
            libfreenect2->initParameters(params[i]);
            libfreenect2->loadP0TablesFromCommandResponse(p0_table, len_p0_table);
            libfreenect2->loadXZTables(&dev.xtable[0], &dev.ztable[0]);
            libfreenect2->loadLookupTable(&dev.lut[0]);
            processors.push_back(libfreenect2);

        }
        else if(strcmp(params[i].pipeline.c_str(),"kde") == 0)
        {
            //initialize kde decoder
            CpuKdeDepthPacketProcessor* kde = new CpuKdeDepthPacketProcessor();
            kde->initParameters(params[1]);
            kde->loadP0TablesFromCommandResponse(p0_table, len_p0_table);
            kde->loadXZTables(&dev.xtable[0], &dev.ztable[0]);
            kde->loadLookupTable(&dev.lut[0]);
            processors.push_back(kde);
        }
        else
        {
            std::cout<<"pipeline "<<params[i].pipeline<<" not implemented\n";
            continue;
        }
        std::string depth_filename = "../data/";
				depth_filename.append(params[i].pipeline);
        std::string conf_filename = "../data/";
				conf_filename.append(params[i].pipeline);
        
        depth_filename.append("_depth_");
        depth_filename.append(params[i].setup_name);
        depth_filename.append("_");
        depth_filename.append(data_set);

        conf_filename.append("_conf_");
        conf_filename.append(params[i].setup_name);        
        conf_filename.append("_");
        conf_filename.append(data_set);

        std::fstream* depth_file = new std::fstream(depth_filename.append(".bin").c_str(), std::ios::out | std::ios::binary | std::ios::app);
        std::fstream* conf_file = new std::fstream(conf_filename.append(".bin").c_str(), std::ios::out | std::ios::binary | std::ios::app);
        depth_files.push_back(depth_file);
        conf_files.push_back(conf_file);
    }

    //decode and save frames to file
    int frame_num;
    int frame_num_offset = 0;
    int num_frames = 10;
    for(frame_num = frame_num_offset; frame_num < num_frames+frame_num_offset; frame_num++)
    {
        std::cout<<"Frame: "<<frame_num<<std::endl;
        unsigned char *buffer;
       
        int len;
        myfile.ReadIrBuffer(&buffer, &len, frame_num);
        for(unsigned int i = 0; i < processors.size(); ++i)
        {
            float *depth, *conf;
            processors[i]->process(buffer, &depth, &conf);
            depth_files[i]->write((char*)depth, 512*424*sizeof(float));
            conf_files[i]->write((char*)conf, 512*424*sizeof(float));
            delete[] depth;
            delete[] conf;
        }

        delete[] buffer;
    }

    //clean up
    for(unsigned int i = 0; i < processors.size(); ++i)
    {
        delete processors[i];
        delete depth_files[i];
        delete conf_files[i];
    }
    std::cout<<"Done!\n";
    return 0;
}
