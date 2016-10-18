#ifndef LIBFREENECT2_DATA_STRUCTURES_H
#define LIBFREENECT2_DATA_STRUCTURES_H

/**
 * Vector class.
 * @tparam ScalarT Type of the elements.
 * @tparam Size Number of elements in the vector.
 */
template<typename ScalarT, int Size>
struct Vec
{
  ScalarT val[Size];
};

/**
 * Matrix class.
 * @tparam ScalarT Eelement type of the matrix.
 */
template<typename ScalarT>
struct Mat
{
private:
  bool owns_buffer; ///< Whether the matrix owns the data buffer (and should dispose it when deleted).
  unsigned char *buffer_; ///< Data buffer of the matrix (row major).
  unsigned char *buffer_end_; ///< End of the buffer (just after the last element).
  int width_;  ///< Number of elements in the matrix.
  int height_; ///< Number of rows in the matrix.
  int x_step;  ///< Number of bytes in one element.
  int y_step;  ///< Number of bytes in one row.

  /**
   * Allocate a buffer.
   * @param width Width of the matrix.
   * @param height Height of the matrix.
   * @param external_buffer If not \c null, use the provided buffer, else make a new one.
   */
  void allocate(int width, int height, unsigned char *external_buffer = 0)
  {
    this->width_ = width;
    this->height_ = height;
    x_step = sizeof(ScalarT);
    y_step = width * x_step;

    owns_buffer = external_buffer == 0;

    if(owns_buffer)
    {
      buffer_ = new unsigned char[y_step * height];
    }
    else
    {
      buffer_ = external_buffer;
    }
    buffer_end_ = buffer_ + (y_step * height);
  }

  void deallocate()
  {
    if(owns_buffer && buffer_ != 0)
    {
      delete[] buffer_;
      owns_buffer = false;
      buffer_ = 0;
      buffer_end_ = 0;
    }
  }

public:
  /** Default constructor. */
  Mat():buffer_(0), buffer_end_(0)
  {
  }

  /**
   * Constructor with locally allocated buffer.
   * @param height Height of the image.
   * @param width Width of the image.
   */
  Mat(int height, int width) : owns_buffer(false), buffer_(0)
  {
    create(height, width);
  }

  /**
   * Constructor with external buffer.
   * @tparam DataT Type of data of the buffer.
   * @param height Height of the image.
   * @param width Width of the image.
   * @param external_buffer Provided buffer.
   */
  template<typename DataT>
  Mat(int height, int width, DataT *external_buffer)
  {
    allocate(width, height, reinterpret_cast<unsigned char *>(external_buffer));
  }

  /** Destructor. */
  ~Mat()
  {
    deallocate();
  }

  /**
   * Get the width of the image.
   * @return Width of the image.
   */
  int width() const
  {
    return width_;
  }

  /**
   * Get the height of the image.
   * @return height of the image.
   */
  int height() const
  {
    return height_;
  }

  /**
   * Construct a new image buffer
   * @param height Height of the new image.
   * @param width Width of the new image.
   */
  void create(int height, int width)
  {
    deallocate();
    allocate(width, height);
  }

  /**
   * Copy image data to the provided matrix.
   * @param other Destination to copy to.
   */
  void copyTo(Mat<ScalarT> &other) const
  {
    other.create(height(), width());
    std::copy(buffer_, buffer_end_, other.buffer_);
  }

  /**
   * Get the image data at the requested point \a x, \a y.
   * @param y Vertical (row) position.
   * @param x Horizontal position.
   * @return Data at the given position.
   */
  const ScalarT &at(int y, int x) const
  {
    return *ptr(y, x);
  }

  /**
   * Get a reference to the image data at the requested point \a x, \a y.
   * @param y Vertical (row) position.
   * @param x Horizontal position.
   * @return Reference to the data at the given position.
   */
  ScalarT &at(int y, int x)
  {
    return *ptr(y, x);
  }

  const ScalarT *ptr(int y, int x) const
  {
    return reinterpret_cast<const ScalarT *>(buffer_ + y_step * y + x_step * x);
  }

  ScalarT *ptr(int y, int x)
  {
    return reinterpret_cast<ScalarT *>(buffer_ + y_step * y + x_step * x);
  }

  /**
   * Get the buffer.
   * @return The buffer.
   */
  unsigned char* buffer()
  {
    return buffer_;
  }

  /**
   * Get the size of the buffer.
   * @return Number of bytes in the buffer.
   */
  int sizeInBytes() const
  {
    return buffer_end_ - buffer_;
  }
};

/**
 * Copy and flip buffer upside-down (upper part to bottom, bottom part to top).
 * @tparam ScalarT Type of the element of the buffer.
 * @param in Source buffer.
 * @param [out] out Destination buffer to be filled with flipped \a in data.
 */
template<typename ScalarT>
void flipHorizontal(const Mat<ScalarT> &in, Mat<ScalarT>& out)
{
  in.copyTo(out);

  typedef unsigned char type;

  int linestep = out.sizeInBytes() / out.height() / sizeof(type);

  type *first_line = reinterpret_cast<type *>(out.buffer()), *last_line = reinterpret_cast<type *>(out.buffer()) + (out.height() - 1) * linestep;


  for(int y = 0; y < out.height() / 2; ++y)
  {
    for(int x = 0; x < linestep; ++x, ++first_line, ++last_line)
    {
      std::swap(*first_line, *last_line);
    }
    last_line -= 2 * linestep;
  }
}

inline int bfi(int width, int offset, int src2, int src3)
{
  int bitmask = (((1 << width)-1) << offset) & 0xffffffff;
  return ((src2 << offset) & bitmask) | (src3 & ~bitmask);
}

  struct Parameters
  {
    float ab_multiplier;
    float ab_multiplier_per_frq[3];
    float ab_output_multiplier;

    float phase_in_rad[3];

    float joint_bilateral_ab_threshold;
    float joint_bilateral_max_edge;
    float joint_bilateral_exp;
    float gaussian_kernel[9];

    float phase_offset;
    float unambigious_dist;
    float individual_ab_threshold;
    float ab_threshold;
    float ab_confidence_slope;
    float ab_confidence_offset;
    float min_dealias_confidence;
    float max_dealias_confidence;

    float edge_ab_avg_min_value;
    float edge_ab_std_dev_threshold;
    float edge_close_delta_threshold;
    float edge_far_delta_threshold;
    float edge_max_delta_threshold;
    float edge_avg_delta_threshold;
    float max_edge_count;

		float kde_sigma_sqr;
		float unwrapping_likelihood_scale;
		float phase_confidence_scale;
		float kde_threshold;
		size_t kde_neigborhood_size;
		size_t num_hyps;

    float min_depth;
    float max_depth;

    DepthPacketProcessor::Parameters::Parameters()
{
  ab_multiplier = 0.6666667f;
  ab_multiplier_per_frq[0] = 1.322581f;
  ab_multiplier_per_frq[1] = 1.0f;
  ab_multiplier_per_frq[2] = 1.612903f;
  ab_output_multiplier = 16.0f;

  phase_in_rad[0] = 0.0f;
  phase_in_rad[1] = 2.094395f;
  phase_in_rad[2] = 4.18879f;

  joint_bilateral_ab_threshold = 3.0f;
  joint_bilateral_max_edge = 2.5f;
  joint_bilateral_exp = 5.0f;

  gaussian_kernel[0] = 0.1069973f;
  gaussian_kernel[1] = 0.1131098f;
  gaussian_kernel[2] = 0.1069973f;
  gaussian_kernel[3] = 0.1131098f;
  gaussian_kernel[4] = 0.1195716f;
  gaussian_kernel[5] = 0.1131098f;
  gaussian_kernel[6] = 0.1069973f;
  gaussian_kernel[7] = 0.1131098f;
  gaussian_kernel[8] = 0.1069973f;

  phase_offset = 0.0f;
  unambigious_dist = 2083.333f;
  individual_ab_threshold  = 3.0f;
  ab_threshold = 10.0f;
  ab_confidence_slope = -0.5330578f;
  ab_confidence_offset = 0.7694894f;
  min_dealias_confidence = 0.3490659f;
  max_dealias_confidence = 0.6108653f;

  edge_ab_avg_min_value = 50.0f;
  edge_ab_std_dev_threshold = 0.05f;
  edge_close_delta_threshold = 50.0f;
  edge_far_delta_threshold = 30.0f;
  edge_max_delta_threshold = 100.0f;
  edge_avg_delta_threshold = 0.0f;
  max_edge_count  = 5.0f;

/*
 * These are parameters for the method described in "Efficient Phase Unwrapping
 * using Kernel Density Estimation", ECCV 2016, Felix Järemo Lawin, Per-Erik Forssen and
 * Hannes Ovren, see http://www.cvl.isy.liu.se/research/datasets/kinect2-dataset/.
 */

    kde_sigma_sqr = 0.0239282226563f; //the scale of the kernel in the KDE, h in eq (13).
    unwrapping_likelihood_scale = 2.0f; //scale parameter for the unwrapping likelihood, s_1^2 in eq (15).
    phase_confidence_scale = 3.0f; //scale parameter for the phase likelihood, s_2^2 in eq (23)
    kde_threshold = 0.5f; //threshold on the KDE output in eq (25), defines the inlier/outlier rate trade-off
    kde_neigborhood_size = 5; //spatial support of the KDE, defines a filter size of (2*kde_neigborhood_size+1 x 2*kde_neigborhood_size+1)
    num_hyps = 2; //number of phase unwrapping hypothesis considered by the KDE in each pixel

  min_depth = 500.0f;
  max_depth = 18750.0f;
}
};

struct P0TablesResponse
{
  uint32_t headersize;
  uint32_t unknown1;
  uint32_t unknown2;
  uint32_t tablesize;
  uint32_t unknown3;
  uint32_t unknown4;
  uint32_t unknown5;
  uint32_t unknown6;

  uint16_t unknown7;
  uint16_t p0table0[512*424]; // row[0] == row[511] == 0x2c9a
  uint16_t unknown8;

  uint16_t unknown9;
  uint16_t p0table1[512*424]; // row[0] == row[511] == 0x08ec
  uint16_t unknownA;

  uint16_t unknownB;
  uint16_t p0table2[512*424]; // row[0] == row[511] == 0x42e8
  uint16_t unknownC;

  uint8_t  unknownD[];
};

