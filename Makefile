TXT = Makefile README.md LICENSE

PYTHON = python
F2PY   = f2py


%.so : %.f90
	$(F2PY) -m $* -c $< --opt='-O3 -fPIC'


all:
	export F2PY=$(F2PY) ; cd core && $(MAKE)

clean:
	rm -f core/*.so core/*/*.so

cleanemacs:
	find . -name '*~' -exec rm {} \;


zip:
	git archive master -onyles_`date '+%d_%m_%Y'`.zip
