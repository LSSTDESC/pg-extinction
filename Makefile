.PHONY: tests

all: tests

clean:
	rm -rf __pycache__

tests:
	env PYTHONPATH=$(CURDIR)/python python3 tests/dustval.py
