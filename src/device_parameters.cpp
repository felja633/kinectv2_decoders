#include "device_parameters.h"
#include "libfreenect2_data_structures.h"
#include <string>

DeviceParametersHandler::DeviceParametersHandler()
{
    ztable = new float[TABLE_SIZE];
    xtable = new float[TABLE_SIZE];
    lut = new float[LUT_SIZE];
}

DeviceParametersHandler::~DeviceParametersHandler()
{
    delete[] ztable;
    delete[] xtable;
    delete[] lut;
}

void DeviceParametersHandler::initializeCameraFromFile(double** dist, double** intr, std::string filename, std::string datastname)
{

	hid_t      s1_tid;                          /* File datatype identifier */
  hid_t      file_id, dataset, space, vlen_tid_dist, vlen_tid_intr;  /* Handles */
	Camera_Device_To_Save cam_struct;
  // Open HDF5 file handle, read only
	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen(file_id, datastname.c_str(),H5P_DEFAULT);

	space = H5Dget_space(dataset);

	s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(Camera_Device_To_Save));
	vlen_tid_dist = H5Tvlen_create(H5T_NATIVE_DOUBLE);
	vlen_tid_intr = H5Tvlen_create(H5T_NATIVE_DOUBLE);

	H5Tinsert(s1_tid, "distortion_coefficients", HOFFSET(Camera_Device_To_Save, distortion_coefficients), vlen_tid_dist);
	H5Tinsert(s1_tid, "num_distortion_coefficients", HOFFSET(Camera_Device_To_Save, num_distortion_coefficients), H5T_NATIVE_INT);
	H5Tinsert(s1_tid, "intrinsic_camera_matrix", HOFFSET(Camera_Device_To_Save, intrinsic_camera_matrix), vlen_tid_intr);
	H5Tinsert(s1_tid, "camera_type", HOFFSET(Camera_Device_To_Save, camera_type), H5T_NATIVE_INT);

	cam_struct.distortion_coefficients.len = 8;
	cam_struct.distortion_coefficients.p = new double[8];

	for(unsigned int i = 0; i < 8; i++)
		((double*)cam_struct.distortion_coefficients.p)[i] = 0.0;

	cam_struct.intrinsic_camera_matrix.len = 9;
	cam_struct.intrinsic_camera_matrix.p = new double[9];
	for(unsigned int i = 0; i < 9; i++)
		((double*)cam_struct.intrinsic_camera_matrix.p)[i] = 0.0;

	H5Dread(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &cam_struct);

	*dist = new double[cam_struct.num_distortion_coefficients];
	double* tmp_dist = (double*)cam_struct.distortion_coefficients.p;
	for(int i = 0; i < cam_struct.num_distortion_coefficients; i++)
		(*dist)[i] = tmp_dist[i];


	*intr = new double[4];
	double* tmp_intr = (double*)cam_struct.intrinsic_camera_matrix.p;
  (*intr)[0] = tmp_intr[0];
	(*intr)[1] = tmp_intr[4];
	(*intr)[2] = tmp_intr[2];
	(*intr)[3] = tmp_intr[5];

	H5Tclose(s1_tid);
  H5Sclose(space);
  H5Dclose(dataset);
  H5Fclose(file_id);

}

void DeviceParametersHandler::initializeRelativePoseFromFile(double** rotation, double** translation, std::string filename)
{
	Relative_Pose_To_Save pose;
	hid_t      s1_tid;                          /* File datatype identifier */
  hid_t      file_id, dataset, space, vlen_tid_rot, vlen_tid_trans;

	file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	dataset = H5Dopen(file_id, "relative_pose",H5P_DEFAULT);

	space = H5Dget_space(dataset);

	s1_tid = H5Tcreate(H5T_COMPOUND, sizeof(Relative_Pose_To_Save));
	vlen_tid_rot = H5Tvlen_create(H5T_NATIVE_DOUBLE);
	vlen_tid_trans = H5Tvlen_create(H5T_NATIVE_DOUBLE);

	H5Tinsert(s1_tid, "rotation", HOFFSET(Relative_Pose_To_Save, rotation), vlen_tid_rot);
	H5Tinsert(s1_tid, "translation", HOFFSET(Relative_Pose_To_Save, translation), vlen_tid_trans);

	pose.rotation.len = 3;
	pose.rotation.p = new double[3];
	pose.translation.len = 3;
	pose.translation.p = new double[3];

  H5Dread(dataset, s1_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, &pose);

	*rotation = new double[3];
	double* tmp_rot = (double*)pose.rotation.p;
	for(unsigned int i = 0; i < 3; i++)
		(*rotation)[i] = tmp_rot[i];

	*translation = new double[3];
	double* tmp_trans = (double*)pose.translation.p;
	for(unsigned int i = 0; i < 3; i++)
		(*translation)[i] = tmp_trans[i];

	H5Tclose(s1_tid);
  H5Sclose(space);
  H5Dclose(dataset);
  H5Fclose(file_id);
}

void DeviceParametersHandler::readCameraParametersFromFile(double** ir_dist, double** rgb_dist, double** ir_intr, double** rgb_intr, double** rotation, double** translation, std::string filename)
{
	initializeCameraFromFile(ir_dist, ir_intr, filename, "leader_calibration_parameters");
	initializeCameraFromFile(rgb_dist, rgb_intr, filename, "follower_calibration_parameters");
	initializeRelativePoseFromFile(rotation, translation, filename);
}

void DeviceParametersHandler::init(std::string filename)
{   
    readCameraParametersFromFile(&ir_distortion_parameters, &rgb_distortion_parameters, &ir_camera_matrix, &rgb_camera_matrix, &rotation, &translation, filename);
    double fx = ir_camera_matrix[0];
    double fy = ir_camera_matrix[1];
    double cx = ir_camera_matrix[2];
    double cy = ir_camera_matrix[3];

    double k1 = ir_distortion_parameters[0];
    double k2 = ir_distortion_parameters[1];
    double p1 = ir_distortion_parameters[2];
    double p2 = ir_distortion_parameters[3];
    double k3 = ir_distortion_parameters[4];

    const double scaling_factor = 8192;
    const double unambigious_dist = 6250.0/3;
    size_t divergence = 0;
    for (size_t i = 0; i < TABLE_SIZE; i++)
    {
      size_t xi = i % 512;
      size_t yi = i / 512;
      double xd = (xi + 0.5 - cx)/fx;
      double yd = (yi + 0.5 - cy)/fy;
      double xu, yu;
      divergence += !undistort(xd, yd, xu, yu, k1, k2, k3, p1, p2);
      xtable[i] = scaling_factor*xu;
      ztable[i] = unambigious_dist/sqrt(xu*xu + yu*yu + 1);
    }

    if (divergence > 0)
      std::cout << divergence << " pixels in x/ztable have incorrect undistortion.\n";

    short y = 0;
    for (int x = 0; x < 1024; x++)
    {
      unsigned inc = 1 << (x/128 - (x>=128));
      lut[x] = y;
      lut[1024 + x] = -y;
      y += inc;
    }
    lut[1024] = 32767;
}

  //x,y: undistorted, normalized coordinates
  //xd,yd: distorted, normalized coordinates
  void DeviceParametersHandler::distort(double x, double y, double &xd, double &yd, double k1, double k2, double k3, double p1, double p2) const
  {
    double x2 = x * x;
    double y2 = y * y;
    double r2 = x2 + y2;
    double xy = x * y;
    double kr = ((k3 * r2 + k2) * r2 + k1) * r2 + 1.0;
    xd = x*kr + p2*(r2 + 2*x2) + 2*p1*xy;
    yd = y*kr + p1*(r2 + 2*y2) + 2*p2*xy;
  }

  //The inverse of distort() using Newton's method
  //Return true if converged correctly
  //This function considers tangential distortion with double precision.
  bool DeviceParametersHandler::undistort(double x, double y, double &xu, double &yu, double k1, double k2, double k3, double p1, double p2) const
  {
    double x0 = x;
    double y0 = y;

    double last_x = x;
    double last_y = y;
    const int max_iterations = 100;
    int iter;
    for (iter = 0; iter < max_iterations; iter++) {
      double x2 = x*x;
      double y2 = y*y;
      double x2y2 = x2 + y2;
      double x2y22 = x2y2*x2y2;
      double x2y23 = x2y2*x2y22;

      //Jacobian matrix
      double Ja = k3*x2y23 + (k2+6*k3*x2)*x2y22 + (k1+4*k2*x2)*x2y2 + 2*k1*x2 + 6*p2*x + 2*p1*y + 1;
      double Jb = 6*k3*x*y*x2y22 + 4*k2*x*y*x2y2 + 2*k1*x*y + 2*p1*x + 2*p2*y;
      double Jc = Jb;
      double Jd = k3*x2y23 + (k2+6*k3*y2)*x2y22 + (k1+4*k2*y2)*x2y2 + 2*k1*y2 + 2*p2*x + 6*p1*y + 1;

      //Inverse Jacobian
      double Jdet = 1/(Ja*Jd - Jb*Jc);
      double a = Jd*Jdet;
      double b = -Jb*Jdet;
      double c = -Jc*Jdet;
      double d = Ja*Jdet;

      double f, g;
      distort(x, y, f, g, k1, k2, k3, p1, p2);
      f -= x0;
      g -= y0;

      x -= a*f + b*g;
      y -= c*f + d*g;
      const double eps = std::numeric_limits<double>::epsilon()*16;
      if (fabs(x - last_x) <= eps && fabs(y - last_y) <= eps)
        break;
      last_x = x;
      last_y = y;
    }
    xu = x;
    yu = y;
    return iter < max_iterations;
  }

