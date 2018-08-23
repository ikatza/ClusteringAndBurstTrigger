CC=g++
CCFLAGS=-Wall `root-config --cflags --libs` -fopenmp -g
EXECS=Module_SNBurst.exe Plotting_SNBurst.exe
# Plotting_SNBurst_Efficiency.exe Module_BetterSNBurst.exe

%.exe: %.C
	$(CC) -o $@ $^ $(CCFLAGS)

all: $(EXECS)

clean:
	rm -rf *.exe

cleanpdf:
	rm -rf *.pdf

