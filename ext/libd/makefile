CC     = g++
CFLAGS = -fPIC -c -O
OBJS   = distances.o invdist.o

libsw.a : $(OBJS)
	ar ru libsw.a $(OBJS)
%.o : %.cpp
	$(CC) $(CFLAGS) -c $<
