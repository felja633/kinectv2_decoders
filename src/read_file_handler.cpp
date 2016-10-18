#include "read_file_handler.h"
#include "hdf5.h"
#include "hdf5_hl.h"
#include <string>
#include <sstream>

using namespace H5;

static std::string Int2String(int val)
{
  std::string res;          // string which will contain the result
  std::ostringstream convert;   // stream used for the conversion
  convert << val;      // insert the textual representation of 'Number' in the characters in the stream
  res = convert.str();
  return res;
}

ReadFileHandler::ReadFileHandler(const std::string in_filename)
{
  //Create new file, delete existing file with the same filename

  mFileId = H5Fopen(in_filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

	//Find number of groups
	hsize_t num_obj;
	H5Gget_num_objs(mFileId, &num_obj);
	number_of_groups = (int)num_obj - 1;
/*
	hsize_t file_size;
	H5Fget_filesize(mFileId, &file_size);*/

}

ReadFileHandler::~ReadFileHandler()
{
  H5Fclose(mFileId);
}


void ReadFileHandler::ReadBuffer(unsigned char** arr, int* length, std::string data_loc)
{
   ReadCharArray(arr, length, mFileId, data_loc);
}

hid_t ReadFileHandler::getGroup(int frame_num)
{
   std::string frame_num_str = Int2String(frame_num);
   std::string frame_group = "/Frame_";
   frame_group.append(frame_num_str);

   hid_t grp_id = H5Gopen2(mFileId, frame_group.c_str(), H5P_DEFAULT);
   return grp_id;
}


void ReadFileHandler::ReadRgbBuffer(unsigned char** arr, int* length, int frame_num)
{
  std::string dataset = "Color";
	int num = frame_num % number_of_groups;
  hid_t group_id = getGroup(num);
  ReadCharArray(arr, length, group_id, dataset);

}

void ReadFileHandler::ReadIrBuffer(unsigned char** arr, int* length, int frame_num)
{

  std::string dataset = "Ir";
	int num = frame_num % number_of_groups;
	//std::cout<<"[ReadFileHandler] read next frame: "<<num<<std::endl;
  hid_t group_id = getGroup(num);

  ReadCharArray(arr, length, group_id, dataset);
  //H5Fclose(group_id);
}

/* Reads a char array to file */
void ReadFileHandler::ReadCharArray(unsigned char** arr, int* length, hid_t identifier, std::string dataset)
{
  //frame_group.append(dataset);
  herr_t      status;
  hsize_t dims[1];

  /* read dataset */
  /* get the dimensions of the dataset */
  H5LTget_dataset_info(identifier, dataset.c_str(), dims, NULL, NULL);
  *length = dims[0];

  *arr = new unsigned char[*length];
   status = H5LTread_dataset_char(identifier, dataset.c_str(), (char*)(*arr));
}


