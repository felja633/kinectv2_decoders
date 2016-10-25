Depth packet processors and plot functions for Kinect v2 phase unwrapping from log files.
The algorithm is described in the paper "Efficient Phase Unwrapping using Kernel
Density Estimation", ECCV 2016, Felix Järemo Lawin, Per-Erik Forssen and 
Hannes Ovren, see http://www.cvl.isy.liu.se/research/datasets/kinect2-dataset/. 

# Setup
```
git clone https://github.com/felja633/kinectv2_decoders.git
```
Download dataset at www.cvl.isy.liu.se/research/datasets/kinect2-dataset/dataset.zip.
Unzip dataset.zip in kinectv2_decoders folder

```
cd kinectv2_decoders
mkdir build
cd build
cmake ..
make
```
# Run
As of now we have three datasets: library, lecture and kitchen. Choose one dataset and run
the code as: 
```
cd kinectv2_decoders/build
./kinectv2_decoders ../parameters/default_parameters.xml dataset
cd ..
python evaluate_decoders.py parameters/default_parameters.xml dataset
```

# Parameters

Paramters are passed in xml-format. At this stage two pipelines are implemented, 
"kde" and "libfreenect2". Each pipeline that is tested should be added in the xml-file.
The user can then add and change the parameters freely.

Example:
```xml
<pipeline name="kde" setup_name="base">
    <Parameters>
        <kde_sigma_sqr>0.0239282226563</kde_sigma_sqr>
        <unwrapping_likelihood_scale>2.0</unwrapping_likelihood_scale>
        <phase_confidence_scale>3.0</phase_confidence_scale>
        <kde_neigborhood_size>5</kde_neigborhood_size>
        <num_hyps>2</num_hyps>
        <min_depth>500.0</min_depth>
        <max_depth>18750.0</max_depth>
    </Parameters>
</pipeline>
```
### Dependencies

The package requires the hdf5 library to parse Kinect v2 log files.

### Mac OSX
Installation instructions:
1. Install HDF5 with brew install (MacOS)
> brew install homebrew/science/hdf5

2. Then add the appropriate paths, e.g.:
export HDF5=/usr/local/Cellar/hdf5/1.8.17 # Replace with your hdf5 installation full path
export PATH=${HDF5}/bin:${PATH}
export DYLD_LIBRARY_PATH=${HDF5}/lib:${DYLD_LIBRARY_PATH}

### Linux
sudo apt-get install libhdf5-dev

# Feedback

For feedback contact Felix Järemo-Lawin <felix.jaremo-lawin@liu.se>
