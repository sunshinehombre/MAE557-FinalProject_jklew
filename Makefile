include Makefile.in
vpath %.f90 $(SRCDIR)

# Targets

all : Makefile Makefile.in
	@make $(BIN)

$(BIN) : $(addprefix $(OBJDIR)/, $(OBJ))
	$(F90) $(F90FLAGS) -o $(BINDIR)/$(BIN) $^

$(OBJDIR)/%.o : %.f90
	$(F90) $(F90FLAGS) -c -o $@ $< $(MODFLAGS)

clean :
	rm -f $(BINDIR)/*
	rm -f $(DATADIR)/[0-9]*
	rm -f $(DATADIR)/*~
	rm -f $(MODDIR)/*.mod
	rm -f $(OBJDIR)/*.o
	rm -f $(SRCDIR)/*~
