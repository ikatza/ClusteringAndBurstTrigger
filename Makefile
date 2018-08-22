CC=g++
CCFLAGS=-Wall `root-config --cflags --libs` -Xpreprocessor -fopenmp -lomp -g
EXECS=Module_SNBurst.exe Plotting_SNBurst.exe MakeSNTheoryDistributions.exe
# Plotting_SNBurst_Efficiency.exe Module_BetterSNBurst.exe
NERROR=100
OPTERROR= -fmax-errors=$(NERROR)

%.exe: %.C
	$(CC) -o $@ $^ $(CCFLAGS) $(OPTERROR)

all: $(EXECS)

clean:
	rm -rf *.exe
	rm -rf *.exe.dSYM

cleanpdf:
	rm -rf *.pdf
