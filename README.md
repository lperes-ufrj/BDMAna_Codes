# BDMAna_Codes
To compile and run the codes

$ cd < directory >

If using GCC, use the switches -DCMAKE_C_COMPILER=$(which gcc) -DCMAKE_CXX_COMPILER=$(which g++)

If using Clang, use the switches -DCMAKE_C_COMPILER=$(which clang) -DCMAKE_CXX_COMPILER=$(which clang++)

$ cmake <swiches as above> -DCMAKE_BUILD_TYPE=RelWithDebInfo .
  
$ make

$ ./BDMAnalysis < input file >
