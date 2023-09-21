SDIR=src
ODIR=obj
IDIR=include
BDIR=bin
DDIR=data

CC=gcc
CFLAGS=-pedantic -Wall -g -I$(IDIR) 
CFLAGS+=-Xpreprocessor -O3
LDFLAGS=-lm
MKDIR_P=mkdir -p
SH=/bin/bash
UNAME_S=$(shell uname -s)

ifeq ($(UNAME_S),Darwin)
	CFLAGS+=-DDARWIN
else
	LDFLAGS+=-lbsd
endif


vpath %.c $(TDIR)
vpath %.c $(SDIR) 

.PHONY: all clean install distclean test cleanOut

TARGET=createDir $(BDIR)/tutorial_digraphs $(BDIR)/tutorial_graphs $(BDIR)/digraphCriticalNodesBF

# Compile C files
$(ODIR)/%.o: %.c
	$(CC) $(CFLAGS) -o $@ -c $<

# Create Directories
createDir:
	${MKDIR_P} ${ODIR}
	${MKDIR_P} ${BDIR}

# Toolkit files
all: $(TARGET) 
	
$(BDIR)/tutorial_digraphs:  $(ODIR)/tutorial_digraphs.o $(ODIR)/graph.o $(ODIR)/utility.o 
	$(CC) $(LDFLAGS) -o $@ $^
	
$(BDIR)/tutorial_graphs:  $(ODIR)/tutorial_graphs.o $(ODIR)/graph.o $(ODIR)/utility.o 
	$(CC) $(LDFLAGS) -o $@ $^
	
$(BDIR)/digraphCriticalNodesBF:  $(ODIR)/digraphCriticalNodesBF.o $(ODIR)/graph.o $(ODIR)/utility.o 
	$(CC) $(LDFLAGS) -o $@ $^
	
# Clean files
clean:
	rm $(ODIR)/*
	rm $(BDIR)/*
	
	
# Clean output
cleanOut:
	rm $(DDIR)/*.DOT
	rm $(DDIR)/*.EDGES
	rm $(DDIR)/*.pdf
	rm $(DDIR)/*.csv
	rm $(DDIR)/*.json
	rm $(DDIR)/*.avgSP

