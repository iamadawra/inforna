
###############################################################
#
#  Modify this part and set the path of your Vienna directory 
#  The Vienna RNA Package can be downloaded at:               
#                       http://www.tbi.univie.ac.at/~ivo/RNA/
#  and has to be configured and installed.

VIENNA          = /home/Vienna_1.6

###############################################################


CXXFLAGS        = -O4 -Wall -g -I$(VIENNA)/include/ViennaRNA 
LDFLAGS         = -L$(VIENNA)/lib -lRNA 
CXX     	= g++

CFLAGS          =

####Files#####

SRCS    = basics.cpp\
          constraints.cpp\
	  struct.cpp\
	  inverse.cpp\
	  energy.cpp\
	  hairpin_energy.cpp\
	  interior_energy.cpp\
	  bulge_energy.cpp\
	  stacking_energy.cpp\
	  multi_energy.cpp\
          end_energy.cpp\
          search.cpp\
          inv_folding_const.cpp

OBJS	= $(SRCS:%.cpp=%.o)

DEPENDFILE	= .depend

EXECUTABLE      = INFO-RNA-2.1.2


### Implicit rules #######

%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) -o $@ $<  
	

###########################

all:	$(DEPENDFILE) $(EXECUTABLE)

include $(DEPENDFILE) 

$(EXECUTABLE):	$(OBJS) 
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LDFLAGS) 


clean:  
	$(RM) $(OBJS) *~ gmon.out $(DEPENDFILE) .gdb_history core

veryclean: clean
	$(RM) $(EXECUTABLE) 

$(DEPENDFILE): 
	(for src in $(SRCS); do $(CXX) $(CXXFLAGS) -MM $${src}; done) > $@
