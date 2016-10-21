#ifndef PROCESSORS_H
#define PROCESSORS_H

#include "libfreenect2_data_structures.h"

class DepthPacketProcessor
{
 public:
  DepthPacketProcessor() {}
  virtual ~DepthPacketProcessor() {}
  virtual void loadP0TablesFromCommandResponse(unsigned char* buffer, size_t buffer_length) = 0;
  virtual void initParameters(Parameters param) = 0;
  virtual void loadXZTables(const float *xtable, const float *ztable) = 0;
  virtual void loadLookupTable(const short *lut) = 0;
  virtual void process(unsigned char* buffer, float** depth_buffer, float** ir_buffer) = 0;

};

class CpuDepthPacketProcessorImpl;

/** Depth packet processor using the CPU. */
class CpuDepthPacketProcessor : public DepthPacketProcessor
{
public:
  CpuDepthPacketProcessor();
  virtual ~CpuDepthPacketProcessor();
  //void setConfiguration(const );

  virtual void loadP0TablesFromCommandResponse(unsigned char* buffer, size_t buffer_length);

  virtual void loadXZTables(const float *xtable, const float *ztable);
  virtual void loadLookupTable(const short *lut);
  virtual void initParameters(Parameters param);
  virtual void process(unsigned char* buffer, float** depth_buffer, float** ir_buffer);
private:
  CpuDepthPacketProcessorImpl *impl_;
};

class CpuKdeDepthPacketProcessorImpl;

/** Depth packet processor using the CPU. */
class CpuKdeDepthPacketProcessor : public DepthPacketProcessor
{
public:
  CpuKdeDepthPacketProcessor();
  virtual ~CpuKdeDepthPacketProcessor();
  //void setConfiguration(const libfreenect2::DepthPacketProcessor::Config &config);

  virtual void loadP0TablesFromCommandResponse(unsigned char* buffer, size_t buffer_length);

  virtual void loadXZTables(const float *xtable, const float *ztable);
  virtual void loadLookupTable(const short *lut);
  virtual void initParameters(Parameters param);
  virtual void process(unsigned char* buffer, float** depth_buffer, float** ir_buffer);
private:
  CpuKdeDepthPacketProcessorImpl *impl_;
};

#endif
