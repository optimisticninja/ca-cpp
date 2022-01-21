#include "debug.h"

#ifdef DEBUG
ostream& dout = cout;
#else
ofstream dev_null("/dev/null");
ostream& dout = dev_null;
#endif
