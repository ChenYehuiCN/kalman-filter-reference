# ESKF example.  The dimensions are compile-time constants of the library,
# so they must be passed to both main.c and ../../eskf/eskf.c.
CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS += -I../../eskf -DESKF_N=3 -DESKF_N_NOM=3 -DESKF_M=6 -DESKF_NW=2
LDLIBS  += -lm
PYTHON  ?= python3

TARGET := eskf_demo
SRCS   := main.c ../../eskf/eskf.c
HDRS   := ../../eskf/eskf.h
CSV    := output/eskf_localization.csv

.PHONY: all run plot demo clean

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS)
	$(CC) $(CPPFLAGS) $(CFLAGS) $(SRCS) -o $@ $(LDLIBS)

run: $(TARGET)
	mkdir -p output
	./$(TARGET) $(CSV)

plot: run
	$(PYTHON) plot.py $(CSV)

demo: plot

clean:
	rm -f $(TARGET)
	rm -rf output
