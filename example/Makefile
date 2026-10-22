# Build/run/plot both examples.
#
#   make            build both C programs
#   make run        run both simulations (writes CSV under */output/)
#   make plot       run both simulations and generate all PNG figures
#   make demo       alias for `make plot`
#   make clean      remove binaries and generated output

.PHONY: all run plot demo clean

all:
	$(MAKE) -C kf_tracking all
	$(MAKE) -C eskf_localization all

run:
	$(MAKE) -C kf_tracking run
	$(MAKE) -C eskf_localization run

plot:
	$(MAKE) -C kf_tracking plot
	$(MAKE) -C eskf_localization plot

demo: plot

clean:
	$(MAKE) -C kf_tracking clean
	$(MAKE) -C eskf_localization clean
