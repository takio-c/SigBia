CC			= g++
LD			= g++
TARGET		= bias.out
CINCLUDES	= -I../matrix/
LINCLUDES	= 
CFLAGS		= -O3
LFLAGS		= -O3

all: $(TARGET)

$(TARGET): bias.o
	$(LD) $(LFLAGS) $(LINCLUDES) -o $(TARGET) bias.o
	./$(TARGET) > data.txt

bias.o: bias.cpp ../matrix/matrix.h
	$(CC) $(CFLAGS) $(CINCLUDES) -c bias.cpp

../matrix/matrix.h: ../matrix/vector.h


clean:
	rm -f $(TARGET) *.o

