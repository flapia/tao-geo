#------------------------------------------------------------------------
#You may need to modify these variables
CC	= gcc #gcc cc
LIBS	= -L$(LIB) -lm -lc
#Options depending on the compiler:
OPTS_linux = -g -w #-Wuninitialized
OPTS_AIX_RS6000 = -g -O3 #-Q -qsrcmsg #-v
OPTS_SUN = -g #-O3 #-xO5 #-O #Sometimes Sun OS has obscure segmentation problems with -O (like those in grav. anom.)
#Chose your compiler:
OPTS	= $(OPTS_linux)
#------------------------------------------------------------------------


