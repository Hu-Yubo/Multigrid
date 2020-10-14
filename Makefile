source = $(wildcard *.cpp)
object = $(patsubst %.cpp, %.o, $(source))

main: $(object)
	g++ -o $@ $^		

.cpp.o:
	g++ -c -o $@ $<

clean :
	-rm -rf $(object)
	-rm -rf main

