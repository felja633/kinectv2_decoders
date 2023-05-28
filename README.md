Depth packet processors and plot functions for Kinect v2 phase unwrapping from log files.

The algorithm is described in the paper *Efficient Phase Unwrapping using Kernel
Density Estimation*, **ECCV 2016**, Felix Järemo Lawin, Per-Erik Forssén and Hannes Ovrén, see http://www.cvl.isy.liu.se/research/datasets/kinect2-dataset/. 

# Setup
```
git clone https://github.com/felja633/kinectv2_decoders.git
```
Download dataset from this [link](https://drive.google.com/file/d/1Jx_hAFySJvfQEWI6PJGpOPVrjKzxg6hK/view?usp=share_link).
Unzip `kinect2_dataset.zip` in the kinectv2_decoders folder

then build the decoders:

```
cd kinectv2_decoders
mkdir build
cd build
cmake ..
make
```

# Run
We provide three datasets: **lecture**, **kitchen**, and **library**. Choose one as `dataset` and run the code as: 
```
cd kinectv2_decoders/build
./kinectv2_decoders ../parameters/default_setup.xml dataset
cd ..
python evaluate_decoders.py test parameters/default_setup.xml dataset
```

To visualize, run code as:
```
python evaluate_decoders.py vis ../parameters/default_setup.xml dataset
```
Toggle between frames using the arrow buttons.

# Parameters
Parameters are passed in xml-format. At this stage two pipelines are implemented, 
**kde** and **libfreenect2**. Each pipeline that is to be tested should be added in the xml-file.
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
# Dependencies

The package requires the **hdf5 library** to parse Kinect v2 log files.

## MacOS X

First Install HDF5 with brew:
~~~~
brew install homebrew/science/hdf5
~~~~
Then add the appropriate paths (in bash):
~~~~
export HDF5=/usr/local/Cellar/hdf5/1.8.17 # Replace with your hdf5 installation full path
export PATH=${HDF5}/bin:${PATH}
export DYLD_LIBRARY_PATH=${HDF5}/lib:${DYLD_LIBRARY_PATH}
~~~~

## Linux
~~~~
sudo apt-get install libhdf5-dev
~~~~
