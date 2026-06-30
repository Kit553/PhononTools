## Workflow





















## Parameters:

INCAR:
  SYSTEM = RELAX
  NWRITE = 2
  PREC = Accurate
  KPAR = 1
  NCORE = 4

  ENCUT = 520
  GGA = PS

  LCHARG = .FALSE.
  LWAVE = .FALSE.

  LASPH = .TRUE.
  ADDGRID = .TRUE.
  EDIFF = 1E-08
  ALGO = Normal
  NELM = 100
  NELMIN = 5
  NELMDL = -5

  EDIFFG = -1E-03
  NSW = 100
  IBRION = 2
  ISIF = 3
  IWAVPR = -5

  ISMEAR = -5
  #SIGMA = 0.04
