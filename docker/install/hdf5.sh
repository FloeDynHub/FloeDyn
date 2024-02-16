# source https://www.wbhatti.org/notes/hdf5-with-cpp-support-debian.html
# wget https://www.hdfgroup.org/downloads/hdf5/source-code/# # not working
cd deps
tar xvfz hdf5-1.14.2.tar.gz
cd hdf5-1.14.2
./configure --prefix=/usr/local --enable-cxx --enable-build-mode=production
make -j 16
make install
echo "/usr/local/lib" > /etc/ld.so.conf.d/h5.conf
ldconfig
