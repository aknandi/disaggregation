all: disaggregation.so

disaggregation.so: disaggregation.cpp
  Rscript --vanilla -e "TMB::compile('disaggregation.cpp', flags = '-Wa,-mbig-obj', safeunload=FALSE)"

clean:
  rm -rf *o