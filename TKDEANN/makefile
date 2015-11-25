#			Linux 
#CC       = g++
#CCOPTS   = -c -I/usr/include/X11R5 -I/usr/X11R6/include -I/usr/include/Motif1.2 -g -DLINUX
#LINK     = g++
#LINKOPTS = -static -L/usr/lib/Motif1.2 -L/usr/lib/X11R5 -lm

#			HP-UX
CC       = g++
CCOPTS   = -c -O
LINK     = g++
LINKOPTS = -lm

UTIL  = utility.o btree.o blk_file.o
RTREE = rtree.o 

.cc.o:
	$(CC) $(CCOPTS) $<

all: gnn

utility.o: utility.cc utility.h 

btree.o: btree.cc btree.h 

blk_file.o: blk_file.cc blk_file.h

est: est.cc $(UTIL) $(RTREE)
	$(LINK) -o est est.cc $(UTIL) $(RTREE) $(LINKOPTS)

plot: plot.cc $(UTIL) $(RTREE)
	$(LINK) -o plot plot.cc $(UTIL) $(RTREE) $(LINKOPTS)

gnn: gnn.cc $(UTIL) $(RTREE)
	$(LINK) -o gnn gnn.cc $(UTIL) $(RTREE) $(LINKOPTS)

topk: topk.cc $(UTIL) 
	$(LINK) -o topk topk.cc $(UTIL) $(LINKOPTS)

basic: basic.cc $(UTIL) $(RTREE)
	$(LINK) -o basic basic.cc $(UTIL) $(RTREE) $(LINKOPTS)

hitiqry: hitiqry.cc $(UTIL) $(RTREE)
	$(LINK) -o hitiqry hitiqry.cc $(UTIL) $(RTREE) $(LINKOPTS)

gendata: gendata.cc $(UTIL) $(RTREE)
	$(LINK) -o gendata gendata.cc $(UTIL) $(RTREE) $(LINKOPTS)

visual: visual.cc utility.o 
	$(LINK) -o visual visual.cc utility.o $(LINKOPTS)

convert: convert.cc utility.o
	$(LINK) -o convert convert.cc utility.o $(LINKOPTS)

rtree.o: rtree.cc rtree.h

bulkload.o: bulkload.cc

main.o: main.cc rtree.h

bulkload: bulkload.o $(UTIL) $(RTREE) 
	$(LINK) -o bulkload bulkload.o $(UTIL) $(RTREE) $(LINKOPTS)

clean:
	rm *.o *stackdump *.exe
