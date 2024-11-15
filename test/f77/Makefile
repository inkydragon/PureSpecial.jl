# -*- mode: makefile -*-

FC:=gfortran
SHARED_LIB_EXT:=.so

FFLAGS+=-I. -shared -fPIC

all: libspecfun$(SHARED_LIB_EXT)

libspecfun$(SHARED_LIB_EXT): specfun.f
	$(FC) $(FFLAGS) -std=legacy -o $@ $^

clean:
	-rm -f libspecfun$(SHARED_LIB_EXT)

# Makefile debugging trick:
# call print-VARIABLE to see the runtime value of any variable
# (hardened against any special characters appearing in the output)
print-%:
	@echo '$*=$(subst ','\'',$(subst $(newline),\n,$($*)))'
