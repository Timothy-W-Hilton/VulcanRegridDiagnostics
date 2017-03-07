SHELL=/bin/sh

PROGRAM=CASA_pergrid_2_perm2

FC=ftn
LD=ftn

# FFLAGS= -O2 -c -C
FFLAGS= -c -C
LDFLAGS= -qopenmp -lpthread -g
APINCL=-I$(PROJ)/local/include
APILIB= -L$(PROJ)/local/lib -lioapi_regrid_tools -lioapi -lnetcdf

CMD=$(PROGRAM).x
SRCDRV=CASA_pergrid_2_perm2
OBJDRV=$(SRCDRV).o

# declare .PHONY so that targets won't be masked by file of same name
.PHONY: all clean clobber debug opt

all:  	$(CMD)

debug: FFLAGS += -DDEBUG -g -O0
debug: $(CMD)

opt: FFLAGS += -O2
opt: $(CMD)

$(CMD):  $(OBJDRV)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJDRV) \
	$(APILIB)
%.o: %.F90
	$(FC) $(FFLAGS) $(APINCL) $<

clean:
	rm -fv  $(OBJDRV) *.o *.mod

clobber: clean
	rm -fv $(CMD)
