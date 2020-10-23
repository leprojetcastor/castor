/*
 +========================================================================+
 |         (c) 2020 - PROPERTY OF ECOLE POLYTECHNIQUE - LGPL 3.0          |
 |________________________________________________________________________|
 |   '&`   |                                                              |
 |    #    |   FILE       : input_output.cpp                              |
 |    #    |   VERSION    : 0.1.0                                         |
 |   _#_   |   AUTHOR(S)  : Matthieu Aussal & Laurent Series              |
 |  ( # )  |   CREATION   : 01.04.2020                                    |
 |  / 0 \  |   LAST MODIF : 31.10.2020                                    |
 | ( === ) |   SYNOPSIS   : Write and read matrix files                   |
 |  `---'  |                                                              |
 +========================================================================+
 */

#include <castor/matrix.hpp>

using namespace castor;

//==============================================================================
int main (int argc, char* argv[])
{
    disp("+======================+");
    disp("|    INPUT - OUTPUT    |");
    disp("+======================+");
    
    matrix<> in = eye(3,4);
    in(2,0) = -M_PI;
    disp("Original matrix :");
    disp(in);
    
    disp("Write and read *.txt file :");
    writetxt("./","matrix_IOtest.txt",in);
    disp(readtxt("./","matrix_IOtest.txt"));
    
    disp("Write and read *.bin file :");
    writebin("./","matrix_IOtest.bin",in);
    disp(readbin("./","matrix_IOtest.bin"));
    
    disp("done !");
    return 0;
}
