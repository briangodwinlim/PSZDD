TDZDD_DIR=.

all: pszdd 

pszdd: pszdd.cpp 
	g++ pszdd.cpp -o pszdd -I$(TDZDD_DIR) -O3

clean:
	rm -f pszdd pszdd.exe *.o