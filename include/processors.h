#ifndef PROCESSORS_H
#define PROCESSORS_H

class CpuDepthPacketProcessorImpl;

/** Depth packet processor using the CPU. */
class CpuDepthPacketProcessor
{
public:
  CpuDepthPacketProcessor();
  ~CpuDepthPacketProcessor();
  void setConfiguration(const libfreenect2::DepthPacketProcessor::Config &config);

  void loadP0TablesFromCommandResponse(unsigned char* buffer, size_t buffer_length);

  void loadXZTables(const float *xtable, const float *ztable);
  void loadLookupTable(const short *lut);

  const char *name() { return "CPU"; }
  void process(const DepthPacket &packet);
private:
  CpuDepthPacketProcessorImpl *impl_;
};

class CpuDepthPacketProcessorImpl;

/** Depth packet processor using the CPU. */
class CpuKdeDepthPacketProcessor
{
public:
  CpuKdeDepthPacketProcessor();
  ~CpuKdeDepthPacketProcessor();
  void setConfiguration(const libfreenect2::DepthPacketProcessor::Config &config);

  void loadP0TablesFromCommandResponse(unsigned char* buffer, size_t buffer_length);

  void loadXZTables(const float *xtable, const float *ztable);
  void loadLookupTable(const short *lut);

  const char *name() { return "CPU"; }
  void process(const DepthPacket &packet);
private:
  CpuKdeDepthPacketProcessorImpl *impl_;
};

#endif
