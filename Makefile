vrlab-diffraction-demo: src/vrlab-diffraction-demo.cpp src/doubleslit-calc.hpp
	@ mkdir -p bin
	g++-11 -Wpedantic -O3 -g -o bin/vrlab-diffraction-demo src/vrlab-diffraction-demo.cpp src/doubleslit-calc.hpp src/lodepng.cpp
	