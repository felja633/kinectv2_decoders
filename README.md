Depth packet processors and plot functions for Kinect v2 phase unwrapping from log files.
The algorithm is described in the paper "Efficient Phase Unwrapping using Kernel
Density Estimation", ECCV 2016, Felix Järemo Lawin, Per-Erik Forssen and 
Hannes Ovren, see http://www.cvl.isy.liu.se/research/datasets/kinect2-dataset/. 

# Compiling

cd kinectv2_decoders
make build
cd build
cmake ..

# Dependencies

The package requires the hdf5 library to parse Kinect v2 log files.

Installation instructions:
1. Install HDF5 with brew install (MacOS)
> brew install homebrew/science/hdf5

2. Then add the appropriate paths, e.g.:
export HDF5=/usr/local/Cellar/hdf5/1.8.17 # Replace with your hdf5 installation full path
export PATH=${HDF5}/bin:${PATH}
export DYLD_LIBRARY_PATH=${HDF5}/lib:${DYLD_LIBRARY_PATH}

# Feedback

For feedback contact Felix Järemo-Lawin <felix.jaremo-lawin@liu.se>
