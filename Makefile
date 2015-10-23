INSTALL_DIR=$(HOME)
SRC = omdtool.cpp param.cpp treader.cpp buffer.cpp

LIB = libomdtool.a
OBJ = $(SRC:.cpp=.o)

CXX = g++
CXXFLAGS = -O -fPIC
ARCHIVE = ar
ARCHFLAG = -rc
LINK = g++
LINKFLAGS =	-O
USRLIB =
SYSLIB =

lib: 	$(OBJ)
	$(ARCHIVE) $(ARFLAGS) $(LIB) $(OBJ)

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $<

install: lib
	@mkdir -p $(INSTALL_DIR)/include/omd
	@mkdir -p $(INSTALL_DIR)/lib
	@cp *.h $(INSTALL_DIR)/include/omd
	@cp *.a $(INSTALL_DIR)/lib


clean:
	-rm -f *.o *~ $(LIB)
