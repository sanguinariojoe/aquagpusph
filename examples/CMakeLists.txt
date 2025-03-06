if(WIN32)
    set(BINARY_NAME AQUAgpusph.exe)
else()
    set(BINARY_NAME AQUAgpusph)
endif()

ADD_SUBDIRECTORY(3D)
ADD_SUBDIRECTORY(2D)
