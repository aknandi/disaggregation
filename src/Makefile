all: disaggregation.so

disaggregation.so: disaggregation.cpp
	"$(R_HOME)/bin/Rscript" --vanilla -e "TMB::compile('disaggregation.cpp')"

clean:
	rm -rf *o