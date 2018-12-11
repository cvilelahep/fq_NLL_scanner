#include "fiTQun.h"

class fq_NLL_scanner: public fiTQun {

 public:
  fq_NLL_scanner(int anpmt = nPMT_max);
  ~fq_NLL_scanner();
  int ScanNLL();
  
};
