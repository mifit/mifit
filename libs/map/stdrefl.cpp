/* Subroutine */ int stdrefl_(long int* ih, long int* ik, long int* il, long int* mult, float* eps, long int* mk, long int* iflg, long int* iflg2, long int* iss,
                              long int* its, long int* nsym) {
  /* System generated locals */
  long int i__1;

  /* Local variables */
  static long int kind, indm, indp, imax, jmax, i, j, k, icent, ispec, jh[3],
                  ip, mk2, ind[3], ikl, isg;

  /*     SUBROUTINE STDREF(IH,IK,IL,MULT,EPS,MK,IFLG,IFLG2) */
  /*     COMMON /SYMTRY/ NSYM,ISS(3,3,24),ITS(3,24) */
  /*  how to use from C: */

  /* #include <mifit/legacy/Xguicryst.h> */

  /* 	int iss[864]; */
  /* 	int its[288]; */
  /* 	int ih,ik,il,mult,mk,iflg,iflg2; */
  /* 	int i,j,k; */
  /* 	float eps; */
  /*       CMapHeaderBase mh; */
  /*       /+ copy its[] and iss[] from rs[] from scan_symmops +/ */
  /* 	for(k=0;k<mh.nsym;k++){ */
  /* 		for(j=0;j<3;j++){ */
  /* 			its[j+k*3]=  ROUND(rs[j][3][k]*24.0); */
  /* 			for(i=0;i<3;i++) */
  /* 				iss[i+3*(j+3*k)]= (int)rs[i][j][k]; */
  /* 		} */
  /* 	} */
  /* 	stdrefl_(&ih,&ik,&il,&mult,&eps,&mk,&iflg, */
  /* 		&iflg2,iss,its,&mh.nsym); */

  /* --- ROUTINE TO FIND STANDARD REFLECTION (MAX VALUE OF PACKED INDICES),
   */
  /* --- MULTIPLICITY, EPSILON AND IDENTIFY REFLECTIONS WITH RESTRICTED */
  /* --- PHASES. UPON EXITING MULT= # TIMES REFLECTION OCCURS IN THE */
  /* --- SPHERE. EPS= MULTIPLE OF THE MEAN INTENSITY (WILSON STATISTICS, */
  /* --- GENERAL REFLECTION), EXPECTED FOR THIS REFLECTION. MK= FLAG */
  /* --- IDENTIFYING RESTRICTED PHASES. IF MK=1, THEN NO RESTRICTIONS. IF */

  /* --- MK > 1, THEN ONE OF THE ALLOWED PHASES IS GIVEN (IN DEGREES) */
  /* --- BY 15*(MK-1), THE OTHER IS 180 DEGREES AWAY.  THE INDICES */
  /*--- RETURNED CORRESPOND TO THE REFLECTION IN ITS STANDARD REPRESENTATION
     .*/
  /* --- ICENT=0 FOR NONCENTROSYMMETRIC SPACE GROUPS, 1 FOR CENTROSYMMETRIC.
   */
  /* --- IFLG= FLAG INDICATING WHETHER FRIEDEL RELATIONSHIP WAS USED */
  /* --- TO PUT THE REFLECTION IN ITS STANDARD REPRESENTATION. IF IFLG=1, */

  /* --- THEN INDICES OF STANDARD REFLECTION ARE RELATED TO THE INPUT */
  /* --- INDICES BY POINT GROUP SYMMETRY ONLY, I.E. THE HAND IS PRESERVED.
   */
  /* --- IF IFLG=-1, FRIEDELS RELATIONSHIP WAS USED. IFLG2=0 FOR ALLOWED */
  /* --- REFLECTIONS, AND IFLG2=1 IF REFLECTION IS FORBIDDEN BY SYMMETRY. */

  /* --- THE ROUTINE IS APPLICABLE TO ALL SPACE GROUPS, INCLUDING ANY */
  /* --- NONSTANDARD SETTINGS. */
  /* --- if mk < 1 then an error occured- i.e.: */
  /* Parameter adjustments */
  its -= 4;
  iss -= 13;

  /* Function Body */
  if (*nsym > 96) {
    *mk = -1;
    return 0;
  }
  ind[0] = *ih;
  ind[1] = *ik;
  ind[2] = *il;
  indp = ((*ih + 128) << 16) + ((*ik + 128) << 8) + *il + 128;
  indm = ((-(*ih) + 128) << 16) + ((-(*ik) + 128) << 8) - *il + 128;
  icent = 0;
  ispec = 0;
  *iflg2 = 0;
  imax = -1;
  *mult = 0;
  *eps = (float)0.;
  *mk = 1;
  if (icent == 1) {
    *mk = 13;
  }
  i__1 = *nsym;
  for (j = 1; j <= i__1; ++j) {
    for (i = 1; i <= 3; ++i) {
      jh[i - 1] = 0;
      for (k = 1; k <= 3; ++k) {
        /* L10: */
        jh[i - 1] += iss[k + (i + j * 3) * 3] * ind[k - 1];
      }
    }
    kind = ((jh[0] + 128) << 16) + ((jh[1] + 128) << 8) + jh[2] + 128;
    /* --- CHECK FOR SYSTEMATIC ABSENCES */
    if (kind == indp) {
      ip = (ind[0] * its[j * 3 + 1] + ind[1] * its[j * 3 + 2] + ind[2] *
            its[j * 3 + 3]) % 24;
      if (ip != 0) {
        goto L90;
      }
    }
    if (kind == indp) {
      *eps += (float)1.;
    }
    if (icent == 1 && kind == indm) {
      *eps += (float)1.;
    }
    if (kind == indp || kind == indm) {
      ++ (*mult);
    }
    /* --- IDENTIFY REFLECTIONS WITH RESTRICTED PHASES */
    if (kind == indm) {
      ispec = j;
    }
    /* --- INSURE THAT FIRST NON-ZERO INDEX IS POSITIVE (HEMISPHERE */
    /* --- WITH H .GE. 0) */
    isg = 1;
    if (jh[0] < 0) {
      goto L40;
    } else if (jh[0] == 0) {
      goto L20;
    } else {
      goto L70;
    }
L20:
    if (jh[1] < 0) {
      goto L50;
    } else if (jh[1] == 0) {
      goto L30;
    } else {
      goto L70;
    }
L30:
    if (jh[2] < 0) {
      goto L60;
    } else if (jh[2] == 0) {
      goto L80;
    } else {
      goto L70;
    }
L40:
    jh[0] = -jh[0];
L50:
    jh[1] = -jh[1];
L60:
    jh[2] = -jh[2];
    isg = -1;
    kind = ((jh[0] + 128) << 16) + ((jh[1] + 128) << 8) + jh[2] + 128;
L70:
    if (kind <= imax) {
      goto L80;
    }
    imax = kind;
    *iflg = isg;
    jmax = j;
L80:
    ;
  }
  *ih = imax / 65536;
  ikl = imax - (*ih << 16);
  *ih += -128;
  *ik = ikl / 256;
  *il = ikl - (*ik << 8) - 128;
  *ik += -128;
  *mult = (*nsym << 1) / *mult;
  if (ispec == 0 || icent == 1) {
    return 0;
  }
  /* --- SPACE GROUP NOT CENTROSYMMETRIC, BUT PHASE IS RESTRICTED, */
  /* --- GET ALLOWED PHASE VALUES. */
  /* --- FOR INPUT (UNTRANSFORMED) INDICES */
  mk2 = (ind[0] * its[ispec * 3 + 1] + ind[1] * its[ispec * 3 + 2] + ind[2]
         * its[ispec * 3 + 3]) / 2;
  mk2 = mk2 % 12 + 1;
  if (mk2 <= 1) {
    mk2 += 12;
  }
  /* --- FOR INDICES IN STANDARD REPRESENTATION */
  *mk = *iflg * (mk2 - 1 - (ind[0] * its[jmax * 3 + 1] + ind[1] * its[jmax *
                                                                      3 + 2] + ind[2] * its[jmax * 3 + 3]));
  *mk = *mk % 12 + 1;
  if (*mk <= 1) {
    *mk += 12;
  }
  return 0;
L90:
  *iflg2 = 1;
  return 0;
} /* stdrefl_ */

