function res = norminvf_as241(p)

  q = p - single(0.5);
  if (abs(q) <= single(0.425))
    r = single(0.180625) - q*q;

    num =         single(5.9109374720e+01);
    num = r*num + single(1.5929113202e+02);
    num = r*num + single(5.0434271938e+01);
    num = r*num + single(3.3871327179e+00);

    den =         single(6.7187563600e+01);
    den = r*den + single(7.8757757664e+01);
    den = r*den + single(1.7895169469e+01);
    den = r*den + single(1.0000000000e+00);

    res = q * num / den;
    return

  else

    if (q < single(0.0))
      r = p;
    else
      r = single(1.0) - p;
    end

    r = sqrt(-log(r));

    if (r <= single(5.0))
      r = r - single(1.6);

      num =         single(1.7023821103e-01);
      num = r*num + single(1.3067284816e+00);
      num = r*num + single(2.7568153900e+00);
      num = r*num + single(1.4234372777e+00);

      den =         single(1.2021132975e-01);
      den = r*den + single(7.3700164250e-01);
      den = r*den + single(1.0000000000e+00);

      res = num / den;

    else
      r = r - single(5.0);

      num =         single(1.7337203997e-02);
      num = r*num + single(4.2868294337e-01);
      num = r*num + single(3.0812263860e+00);
      num = r*num + single(6.6579051150e+00);

      den =         single(1.2258202635e-02);
      den = r*den + single(2.4197894225e-01);
      den = r*den + single(1.0000000000e+00);

      res = num / den;
    end

    if (q < single(0.0))
      res = - res;
    end

    return
  end
end