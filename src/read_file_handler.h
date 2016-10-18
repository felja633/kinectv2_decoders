#include "H5Cpp.h"
#include <iostream>
#include <string>
#include <vector>

#ifndef READ_FILE_HANDLER_H
#define READ_FILE_HANDLER_H

class ReadFileHandler {

public:
  ReadFileHandler(const std::string filename);
  ~ReadFileHandler();

  void ReadRgbBuffer(unsigned char** arr, int* length, int frame);
  void ReadIrBuffer(unsigned char** arr, int* length, int frame);
  void ReadBuffer(unsigned char** arr, int* length, std::string data_loc);  
  
	int number_of_groups;
private:
  hid_t getGroup(int frame_num);
  void ReadCharArray(unsigned char** arr, int* length, hid_t identifier, std::string dataset);
  hid_t mFileId;
};



#endif

