all: disaggregation.so

disaggregation.so: disaggregation.cpp
  Rscript --vanilla -e "TMB::compile('disaggregation.cpp')"

clean:
  rm -rf *o