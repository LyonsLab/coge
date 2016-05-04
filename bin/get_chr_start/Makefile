app=get_chr_start
CXX=g++

SRCS=$(app).cpp chr_index.cpp
OBJS=$(subst .cpp,.o,$(SRCS))

all: $(app)

$(app): $(OBJS)
	$(CXX) -o $(app) $(OBJS)

$(app).o: $(app).cpp chr_index.hpp

chr_index.o: chr_index.cpp chr_index.hpp

debug: 
	g++ -g -o $(app) $(app).cpp

clean:
	rm $(app) $(OBJS)
