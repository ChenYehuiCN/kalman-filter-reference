# Linear KF example.  The dimensions are compile-time constants of the
# library, so they must be passed to both main.c and ../../kf/kf.c.
CC      ?= cc
CFLAGS  ?= -O2 -std=c11 -Wall -Wextra
CPPFLAGS += -I../../kf -DKF_N=6 -DKF_M=2 -DKF_C=2
LDLIBS  += -lm
PYTHON  ?= python3

TARGET := kf_demo
SRCS   := main.c ../../kf/kf.c
HDRS   := ../../kf/kf.h
CSV    := output/kf_tracking.csv

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
