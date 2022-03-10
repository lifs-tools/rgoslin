CC = g++
opt = -O3 -std=c++11 -Wall
dopt = -O0 -std=c++11 -Wall -g -rdynamic
inc = -I /usr/include -I ../.
lib_dir = -L /usr/lib -L ../.
libs = -l cppGoslin

main: lipid_name_parser.cpp simple_statistics.cpp
	${CC} ${opt} ${inc} -o lipid_name_parser lipid_name_parser.cpp ${lib_dir} ${libs}
	${CC} ${opt} ${inc} -o simple_statistics simple_statistics.cpp ${lib_dir} ${libs}
	${CC} ${opt} ${inc} -o swiss_lipids_statistics swiss_lipids_statistics.cpp ${lib_dir} ${libs}

clean:
	rm -f lipid_name_parser
	rm -f simple_statistics
	rm -f swiss_lipids_statistics
