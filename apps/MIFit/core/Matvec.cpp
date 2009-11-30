#define PRODUCT(a, b) ((a)*(b))
void matvec(double *m, double *v, double *rv)
{
    register double *r;
    double work[3]; /* workarea for forming inner product sums */
    register int i;
    register double tally;

    r = &work[0];

    for (i = 0; i < 3; i++)
    {
        tally = PRODUCT(*m++, *v++) ;
        tally += PRODUCT(*m++, *v++);
        *r++ = tally + PRODUCT(*m++, *v++);
        v -= 3;
    }

    /* copy into result array */
    r = &work[0];
    *rv++ = *r++;
    *rv++ = *r++;
    *rv   = *r;
}

