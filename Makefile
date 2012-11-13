CC			= g++
LD			= g++
TARGET		= bias.out
CINCLUDES	= -I./matrix/
LINCLUDES	= 
CFLAGS		= -O3
LFLAGS		= -O3

all: $(TARGET)

$(TARGET): bias.o
	$(LD) $(LFLAGS) $(LINCLUDES) -o $(TARGET) bias.o

bias.o: bias.cpp
	$(CC) $(CFLAGS) $(CINCLUDES) -c bias.cpp

clean:
	rm -f $(TARGET) *.o

