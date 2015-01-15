- generate Makefile & compile & run:

    cmake .
    make
    ./hello
    
- configure & generate Makefile & compile & run:

    ccmake .
    [configure settings, press 'c', press 'g']
    make
    ./hello
    
- select a specific c++ compiler (must be installed, otherwise cmake fails!) & generate Makefile & compile & run:

    [set environment variable CXX to desired c++ compiler - using bash, "export CXX=icpc" will select the intel compiler]
    cmake .
    make
    ./hello