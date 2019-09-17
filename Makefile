TXT = Makefile README.md LICENSE

PYTHON = python
F2PY   = f2py


%.so : %.f90 ; $(F2PY) -m $* -c $< --opt='-O3 -fPIC'


all:
	export F2PY=$(F2PY) ; cd core && $(MAKE)

clean:
	rm -f core/*.so core/*/*.so


zip:	
	cd .. ; zip -r Nyles/nyles_`date '+%d_%m_%Y'`.zip Nyles/Makefile Nyles/README* Nyles/activate.* Nyles/LICENSE Nyles/core  -x \*.so \*~ \*.pyc \*.nc \.* \*checkpoint* \*pycache* \*.rst \*.png


