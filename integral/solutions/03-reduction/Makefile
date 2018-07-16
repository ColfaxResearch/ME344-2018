CC = icpc
LD = icpc
CFLAGS= -fimf-precision=low -O2 -g -xHost -qopt-zmm-usage=high -qopt-report=5 -qopenmp -par-affinity:compact,1,granularity=fine
LDFLAGS=-qopenmp

OBJECTS = src/integral.o src/main.o

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
	$(info Copying code from solutions/01-critical/)
	cp -v solutions/01-critical/Makefile ./
	cp -v solutions/01-critical/src/* ./src/

solution-02:
	$(info Copying code from solutions/02-atomic/)
	cp -v solutions/02-atomic/Makefile ./
	cp -v solutions/02-atomic/src/* ./src/

solution-03:
	$(info Copying code from solutions/03-reduction/)
	cp -v solutions/03-reduction/Makefile ./
	cp -v solutions/03-reduction/src/* ./src/

solution-04:
	$(info Copying code from solutions/04-private/)
	cp -v solutions/04-private/Makefile ./
	cp -v solutions/04-private/src/* ./src/

clean: 
	rm -f $(OBJECTS) $(TARGET) src/*.optrpt *~ src/*~

