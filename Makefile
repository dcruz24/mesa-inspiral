PYTHON ?= python

.PHONY: test unittest lint check

test:
	$(PYTHON) -m pytest

unittest:
	$(PYTHON) -m unittest discover -s tests -v

lint:
	$(PYTHON) -m pylint main.py src tests

check: test lint
