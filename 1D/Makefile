source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))

main: $(object)
	g++ -g -o $@ $^ 		

.cpp.o:
	g++ -g -c -o $@ $<

clean :
	-rm -rf $(object)
	-rm -rf main

