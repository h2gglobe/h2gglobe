#ifndef _utils_h
#define _utils_h

#include <iostream>
#include <cstdlib>
/** running interactively with ROOT, the standard assert creates a confusing situation... */
#define ASSERT(condition) \
  if (! (condition)) { \
    std::cerr << "Assertion '" << #condition << "' failed in " << __FILE__ << ":" << __LINE__ << std::endl; \
    exit(1); \
  }

#endif
