#ifndef logn_HH
#define logn_HH

#include <cmath>

double logn(double x, unsigned p=2) {
  static const double logx0(-100.0);
  static const double x0(exp(logx0));
  static const double x02(x0*x0);

  if(x>=x0) return log(x);

  if(p==0) return logx0;
  if(p==1) return x/x0-1.0+logx0;
  return -0.5*x*x/x02+2.0*x/x0-1.5+logx0;
}

#endif
