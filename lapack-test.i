local _timer_stamp;
func timer_start {
  extern _timer_stamp;
  _timer_stamp = array(double, 3);
  timer, _timer_stamp;
}
func timer_elapsed(count) {
  extern _timer_stamp;
  elapsed = _timer_stamp;
  timer, elapsed;
  elapsed -= _timer_stamp;
  if (! is_void(count)) elapsed /= count;
  if (am_subroutine()) {
    cpu = elapsed(1);
    sys = elapsed(2);
    all = elapsed(3);
    write, format="cpu=%g (%.0f%%), system=%g, wall=%g\n",
      cpu, 1e2*(all > 0.0 ? cpu/all : 1.0), sys, all;
  } else {
    return elapsed;
  }
}

local _gettimeofday;
func gettimeofday(nil)
{
  if (is_void(_gettimeofday)) {
    require, "dlwrap.i";
    dll = dlopen();
    _gettimeofday = dlwrap(dll,DL_INT,"gettimeofday",DL_LONG_ARRAY,DL_ADDRESS);
  }
  result = array(long, 2);
  if (_gettimeofday(result, 0) != 0) {
    error, "error in gettimeofday";
  }
  return result;
}

//plug_dir,".";
//include,"lapack.i";
//func lpk_gesv(..) {return 0.0;}

func lpk_test_warmup
{
  random_n, 1000*1000;
}

func lpk_test_dims(a,b)
{
  return (numberof(a) == numberof(b) && allof(a == b));
}

func lpk_test_rand(type, dims, vmin, vmax)
{
  if (is_void(vmin)) vmin = 0.0;
  if (is_void(vmax)) vmax = vmin + 1.0;
  x = (vmax - vmin)*random(dims) + vmin;
  if (type == double) return x;
  if (type == complex) return x + ((vmax - vmin)*random(dims) + vmin)*1i;
  return type(x);
}

func lpk_test_blas1(size=, nloops=)
{
  type = float;
  dims = 30;
  x0 = lpk_test_rand(type, dims, -0.5, 0.5);
  y0 = lpk_test_rand(type, dims, -0.5, 0.5);

  alpha = [0, -1, 1, 1.5, 2, -2, 0.01];

  for (k = 1; k <= numberof(alpha); ++k) {
    x = x0;
    lpk_scal, alpha(k), x;
    write, format="lpk_scal, %g, x; ---> %g\n",
      alpha(k), max(abs(x - alpha(k)*x0));
  }

  for (k = 1; k <= numberof(alpha); ++k) {
    x = x0;
    lpk_scal, alpha(k), x, ::-1;
    write, format="lpk_scal, %g, x, ::-1; ---> %g\n",
      alpha(k), max(abs(x - alpha(k)*x0));
  }

}

func lpk_test_dot(size=, nloops=)
{
  if (is_void(size)) size = 1000;
  if (is_void(nloops)) nloops = 100;
  x = random(size);
  y = random(size);

  if (nloops >= 1) {
    lpk_test_warmup;

    timer_start;
    for (k=1;k<=nloops;++k) {
      res1 = sum(x*y);
    }
    write, format="%s", "results for sum(x*y): ";
    timer_elapsed, nloops;

    timer_start;
    for (k=1;k<=nloops;++k) {
      res2 = lpk_dot(x, y);
    }
    write, format="%s", "results for lpk_dot:  ";
    timer_elapsed, nloops;
  }

  if (nloops < 1) {
    res1 = sum(x*y);
    res2 = lpk_dot(x, y);
  }

  write, format = "difference (ratio): %g (%g)\n",
    abs(res1 - res2), abs(res1 - res2)/res1;
}

func lpk_test_gemm(size=, nloops=)
{
  n = 6400;
  a = random(n,n);

  b = random(n,3);
  c = array(double, n, 3);
  lpk_gemm, LPK_NO_TRANS, LPK_NO_TRANS, 1.0, a, b, 0.0, c;
  cp = a(,+)*b(+,);
  stat, c-cp;
    lpk_gemm, LPK_TRANS, LPK_NO_TRANS, 1.0, a, b, 0.0, c;
  cp = a(+,)*b(+,);
  stat, c-cp;

  b = random(3,n);
  lpk_gemm, LPK_NO_TRANS, LPK_TRANS, 1.0, a, b, 0.0, c;
  cp = a(,+)*b(,+);
  stat, c-cp;
  lpk_gemm, LPK_TRANS, LPK_TRANS, 1.0, a, b, 0.0, c;
  cp = a(+,)*b(,+);
  stat, c-cp;
}

LPK_BENCH_SIZE = [50, 100, 200, 300, 500, 1000, 2000, 3000, 5000, 10000, 20000];
func lpk_bench_gesv(size=, ntrials=, nsecs=)
{
  if (is_void(size)) size = LPK_BENCH_SIZE;
  if (is_void(ntrials)) ntrials = 5; // minimum number of trials
  if (is_void(nsecs)) nsecs = 10.0; // minimum number of seconds

  hline = "---------------------------------------------------------------------";

  for (k = 1; k <= numberof(size); ++k) {
    if (k == 1) {
      write, format = "%s\n%s\n%s\n%s\n",
        hline,
        "                 Yorick                    Lapack",
        "  size           LUsolve                  lpk_gesv",
        hline;
    }
    n = size(k);
    a = random(n,n) - 0.5;
    x = random(n);
    b = a(,+)*x(+); // FIXME: use GEMV

    elapsed = array(double, 3);
    for (pass = 1; pass <= 2; ++pass) {
      op = (pass == 1 ? LUsolve : lpk_gesv);
      s1 = s2 = array(double, 3);
      cnt = 0;
      do {
        timer, elapsed;
        start = elapsed;
        xp = op(a, b);
        timer, elapsed;
        t = elapsed - start;
        s1 += t;
        s2 += t*t;
      } while  (++cnt < ntrials || s1(3) < nsecs);
       // FIXME: take the median
      mean = s1/cnt;
      stdev = sqrt((s2 - s1*s1/cnt)/(cnt - 1));
      if (pass == 1) {
        write, format = "%6d  %.2E +/- %.0E (%3.0f%%) ",
          n, mean(3), stdev(3), 100.0*mean(1)/mean(3);
        ref = mean(3);
      } else {
        write, format = "  %.2E +/- %.0E (%3.0f%%) [%2.1f]\n",
          mean(3), stdev(3), 100.0*mean(1)/mean(3), ref/mean(3);
      }
    }
  }
  write, format = "%s\n", hline;
}

func lpk_test_gesv(n, nloops=)
{
  /* test of multi-dimension features */
  a = random(3,4,5,3,4,5)  - 0.5;
  x = random(3,4,5,6,7) - 0.5;
  b = (a(,,,*))(..,+)*(x(*,,))(+,..);
  x2 = lpk_gesv(a, b);
  max(abs(x2 - x));

  if (is_void(n)) n = 50;
  if (is_void(nloops)) nloops = 1;
  a = random(n,n) - 0.5;
  x = random(n);
  timer_start;
  b = a(,+)*x(+);
  write, format="%s", "Yorick matrix multiplication: ";
  timer_elapsed, 1;
  a0 = a;
  b0 = b;

  if (nloops >= 1) {
    lpk_test_warmup;

    //t0 = gettimeofday();
    timer_start;
    for (k=1;k<=nloops;++k) {
      x1 = LUsolve(a, b);
    }
    write, format="%s", "results for LUsolve:  ";
    timer_elapsed, nloops;
    //t1 = gettimeofday();
    //write, format="total time: %g s/iterations\n",
    //  sum((t1 - t0)*[1.0, 1e-6])/nloops;

    //t0 = gettimeofday();
    timer_start;
    for (k=1;k<=nloops;++k) {
      x2 = lpk_gesv(a, b);
    }
    write, format="%s", "results for lpk_gesv: ";
    timer_elapsed, nloops;
    //t1 = gettimeofday();
    //write, format="total time: %g s/iterations\n",
    //  sum((t1 - t0)*[1.0, 1e-6])/nloops;
  }

  if (nloops < 1) {
    x1 = LUsolve(a, b);
    x2 = lpk_gesv(a, b);
  }
  max(abs(a - a0));
  max(abs(b - b0));
  max(abs(x1 - x));
  max(abs(x2 - x));
  max(abs(x2 - x1));
}

func lpk_test0(m)
{
  if (is_void(m)) m = 3000;
  a = random(2*m,m) - 0.5;
  a = a(+,)*a(+,);
  for (k=1;k<=100;++k) {
    b = a;
    lpk_potrf, LPK_UPPER, b;
  }

}


//include,"linalg.i";

func lpk_test1(from=, to=, step=, type=)
{
  if (is_void(from)) from = 1;
  if (is_void(to)) to = 200;
  if (is_void(step)) step = 1;
  if (is_void(type)) type = float;

  swrite, format="From : %3d  To : %3d Step = %3d\n", from, to, step;

  for (m = from; m <= to; m += step) {

    swrite, format="M = %6d : ", m;
    b = array(type, m, m);

    for (pass = 1; pass <= 2; ++pass) {
      a = array(type, m, m);
      a(1:m*m:m+1) = 8.0 + random(m);
      if (pass == 1) {
        trans = LPK_TRANS;
        uplo = LPK_UPPER;
        alpha = 1.0;
        beta = 0.0;
        if (type == complex) {
          for (j = 2; j <= m; ++j) {
            a(1:j-1,j).re = random(j-1) - 0.5;
            a(1:j-1,j).im = random(j-1) - 0.5;
          }
        } else {
          for (j = 2; j <= m; ++j) {
            a(1:j-1,j) = random(j-1) - 0.5;
          }
        }
      } else {
        trans = LPK_NO_TRANS;
        uplo = LPK_LOWER;
        alpha = 0.0;
        beta = 0.0;
        if (type == complex) {
          for (j = 1; j < m; ++j) {
            a(j+1:m,j).re = random(m-j) - 0.5;
            a(j+1:m,j).im = random(m-j) - 0.5;
          }
        } else {
          for (j = 1; j < m; ++j) {
            a(j+1:m,j) = random(m-j) - 0.5;
          }
        }
      }
      pm, a, format="%6.3f";
      lpk_syrk, uplo, trans, alpha, a, beta, b;
      pm, b, format="%6.3f";
      lpk_potrf, uplo, b;
      pm, b, format="%6.3f";
    }
  }
}

func lpk_test_gesvd(m,n,type=,full=,dac=)
{
  local s, u, vt, s0, u0, v0t;

  if (is_void(type)) type = double;
  if (type == complex) {
    a = random_n(m,n) + 1i*random_n(m,n);
  } else if (type == double) {
    a = random_n(m,n);
  } else if (type == float) {
    a = float(random_n(m,n));
  } else {
    error, "TYPE must be: float, double or complex";
  }

  a0 = a;

  t = (dac ? lpk_gesdd : lpk_gesvd)(a,s,u,vt);
  if (t !=0) error;

  //write, format="checking %s: ", "A";
  //stat, a - a0;

  if (type != complex) {
    s0 = SVdec(a0, u0, v0t);
    write, format="checking %s: ", "S";
    stat, s - s0;
  }

  p = min(m,n);
  t = u(,1:p)(,+)*(s*vt(1:p,))(+,);

  write, format="checking %s: ", "A";
  stat, a - t;


  //write, format="checking %s: ", "U";
  //stat, u - u0;
  //
  //write, format="checking %s: ", "V";
  //stat, vt - v0t;

}
