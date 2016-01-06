all: task01 task02 task03 task04 task05 task06

task01: task0102.cpp
	g++ task0102.cpp -DTASK=1 -o task01

task02: task0102.cpp
	g++ task0102.cpp -DTASK=2 -o task02

task03: task03.cpp
	g++ task03.cpp -o task03

task04: task04.cpp
	mpic++ --std=c++11 task04.cpp -o task04

task05: task05.cpp
	mpic++ --std=c++11 task05.cpp -o task05

task06: task06.cpp
	mpic++ --std=c++11 task06.cpp -o task06

clean:
	rm task01 task02 task03 task04 task05 task06
