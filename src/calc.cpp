double unifRand() {
  return rand() / ( (double) RAND_MAX + 1);
}

double log_exp_x_plus_exp_y(double x, double y) {

    double result;
    if ( x == log(0) && y==log(0)) result = log(0);
    else if (x== -log(0) && y==-log(0)) result = -log(0);
    else if ( x - y >= 100 ) result = x;
    else if ( x - y <= -100 ) result = y;
    else {
      if (x > y) {
      result = y + log( 1 + exp(x-y) );
      }
      else result = x + log( 1 + exp(y-x) );
    }
    return result;
}       
