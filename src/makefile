cc = g++
object = main.o comp.o decomp.o tools.o setting.o
main:
	$(cc) -o main main.cpp -pthread
ifeq ("x","y")
main.o: main.cpp
	$(cc) -c main.cpp
comp.o: comp.cpp
	$(cc) -c comp.cpp
decomp.o: decomp.cpp
	$(cc) -c decomp.cpp
tools.o: tools.cpp
	$(cc) -c tools.cpp
setting.o: setting.cpp
	$(cc) -c setting.cpp
endif
clean:
	rm main
