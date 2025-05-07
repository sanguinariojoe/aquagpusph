# CMake generated Testfile for 
# Source directory: /home/yanez/src/c++/branch_yanez/tests/ExternalTool
# Build directory: /home/yanez/src/c++/branch_yanez/build/tests/ExternalTool
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(ExternalTool "/usr/bin/bash" "/home/yanez/src/c++/branch_yanez/build/tests/ExternalTool/run.sh")
set_tests_properties(ExternalTool PROPERTIES  DEPENDS "exttooltest" _BACKTRACE_TRIPLES "/home/yanez/src/c++/branch_yanez/cMake/Test.cmake;18;add_test;/home/yanez/src/c++/branch_yanez/tests/ExternalTool/CMakeLists.txt;18;aquagpusph_conf_test;/home/yanez/src/c++/branch_yanez/tests/ExternalTool/CMakeLists.txt;0;")
