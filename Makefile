FFLAGS = -O3 -fbounds-check
FC = gfortran
SOURCE = \
fragment.f block.f common.f cputim.f event.f fpoly.f fragmt.f histog.f \
infrag.f input.f intgrt.f merge.f output.f ran2.f search.f sqrtf.f

OBJECTS = $(SOURCE:.f=.o)

fragment:	$(OBJECTS)
	$(FC) $(FCLAGS) $(OBJECTS) -o fragment

print:
	@- \rm -f FRAGMENT.TEXT
	@cat $(SOURCE) > FRAGMENT.TEXT

clean:
	rm -f $(OBJECTS) fragment
