CC = icpc
LD = icpc
CFLAGS=-O2 -g -xHost -qopt-zmm-usage=high -qopt-report=5 -qopenmp -par-affinity=compact,1,granularity=fine
LDFLAGS=-qopenmp -mkl

OBJECTS = src/lu.o src/main.o

TARGET=app

all: $(OBJECTS)
	$(LD) -o $(TARGET) $(OBJECTS) $(LDFLAGS)

run:
	./$(TARGET)

solution-00:
	$(info Copying code from solutions/00-initial/)
	cp -v solutions/00-initial/Makefile ./
	cp -v solutions/00-initial/src/* ./src/

solution-01:
	$(info Copying code from solutions/01-permuted/)
	cp -v solutions/01-permuted/Makefile ./
	cp -v solutions/01-permuted/src/* ./src/

solution-02:
	$(info Copying code from solutions/02-parallel-for/)
	cp -v solutions/02-parallel-for/Makefile ./
	cp -v solutions/02-parallel-for/src/* ./src/

solution-03:
	$(info Copying code from solutions/03-tasks/)
	cp -v solutions/03-tasks/Makefile ./
	cp -v solutions/03-tasks/src/* ./src/

solution-04:
	$(info Copying code from solutions/04-mkl/)
	cp -v solutions/04-mkl/Makefile ./
	cp -v solutions/04-mkl/src/* ./src/

clean: 
	rm -f $(OBJECTS) $(TARGET) src/*.optrpt *~ src/*~

