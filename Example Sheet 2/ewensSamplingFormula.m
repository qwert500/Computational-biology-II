function ewensSamplingFormula=ewensSamplingFormula(theta,n)
for i=1:n
  if i==1
     ewensSamplingFormula=1;
  else
  ewensSamplingFormula=ewensSamplingFormula+theta/(i-1+theta);
  end
end
  