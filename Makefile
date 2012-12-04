CC			= g++
LD			= g++
TARGET		= bias.out
SRCS		= bias.cpp
OBJS		= $(SRCS:%.cpp=%.o)
DEPS		= $(SRCS:%.cpp=%.d)
CINCLUDES	= -I../matrix/
LINCLUDES	= 
CFLAGS		= -g
LFLAGS		= -O3

all: $(TARGET)

-include $(DEPS)

$(TARGET): $(OBJS)
	$(LD) $(LFLAGS) $(LINCLUDES) -o $@ $^
	./$(TARGET) > data.txt

%.o: %.cpp
	$(CC) $(CFLAGS) $(CINCLUDES) -c -MMD -MP $<

clean:
	rm -f $(TARGET) $(OBJS) $(DEPS)

