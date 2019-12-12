all: minSurf

DCO_INC=../dco_cpp/include
DCO_LIB=../dco_cpp/lib

minSurf:
	g++ -O3 -g -I${DCO_INC} -L${DCO_LIB} main.cpp -ldcoc -o minSurf.out

run: minSurf
	export NAG_KUSARI_FILE=../nag_key.txt && ./minSurf.out | tee z

clean:
	rm -f minSurf.out z
