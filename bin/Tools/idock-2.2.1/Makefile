CC=clang++ -O2

bin/idock: obj/io_service_pool.o obj/safe_counter.o obj/array.o obj/atom.o obj/scoring_function.o obj/receptor.o obj/ligand.o obj/result.o obj/random_forest.o obj/random_forest_x.o obj/random_forest_y.o obj/main.o
	${CC} -o $@ $^ -pthread -L${BOOST_ROOT}/lib -lboost_system -lboost_program_options -lboost_filesystem -lboost_thread

obj/%.o: src/%.cpp
	${CC} -o $@ $< -c -std=c++11 -Wall -Wno-unused-local-typedef -I${BOOST_ROOT}

clean:
	rm -f bin/idock obj/*.o
