cmake . \
    -DAQUAGPUSPH_3D=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DOPENCL_ROOT_DIR=/etc/OpenCL \
    -DPYTHON_LIBRARIES=$LIBRARY_PATH \
    -DPYTHON_INCLUDE_DIRS=$INCLUDE_PATH \
    -DPYTHON_LIBRARY=$LIBRARY_PATH/libpython$PY_VER.so \
    -DPYTHON_INCLUDE_DIR=$INCLUDE_PATH/python$PY_VER \
    -DCMAKE_INSTALL_PREFIX=$PREFIX \
    -DCMAKE_INSTALL_RPATH=$LIBRARY_PATH \
    -DCMAKE_PREFIX_PATH=$PREFIX \
    -DAQUAGPUSPH_USE_VTK:BOOL=ON \
    -DAQUAGPUSPH_USE_NCURSES:BOOL=OFF \
    -DLLVM_LIBRARY=$LIBRARY_PATH/libLLVM-3.3.so

# HACK to fix the link rpath directory
sed -e "s;conda-bld/work/lib;envs/_build/lib;" -i $SRC_DIR/src/CMakeFiles/AQUAgpusph.dir/link.txt
make
make install
