STANDART = -std=c++11
FLAGS = -Wall -Werror -fopenmp
NAME = prod

all: prod.cpp
	g++ ${FLAGS} ${STANDART} prod.cpp -o ${NAME}	
clean:
	rm -f bctopo.msh
	rm -f coordinate.msh
	rm -f mesh.txt
	rm -f topo.msh
	rm -f ${NAME}
	rm -f decision.txt
