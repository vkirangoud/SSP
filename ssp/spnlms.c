/* SPNLMS     11/13/85 */
/* IMPLEMENTS NLMS ALGORITHM B(K+1)=B(K)+2*MU*E*X(K)/((L+1)*SIG) */
/* X(O:N-1)=DATA VECTOR      INPUT SENT     OUTPUT RETURNED */
/* D(0:N-1)=DESIRED SIGNAL VECTOR */
/* N SPECIFIES NUMBER OF DATA POINTS IN X AND D */
/* B(0:L)=ADAPTIVE COEFFICIENTS OF LTH ORDER FIR FILTER */
/* MU=CONVERGENCE PARAMETER - DECLARE REAL */
/* SIG=INPUT SIGNAL POWER ESTIMATE - UPDATED INTERNALLY */
/* AL=FORGETTING FACTOR   SIG(K)=AL*(X(K)**2)+(1-AL)*SIG(K-1) */
/* IERROR=0      NO ERRORS DETECTED */
/*        1      INVALID ORDER   L<0 */
/*        2      INVALID CONVERGENCE PARAMETER   MU<=0 OR >=1 */
/*        3      INPUT POWER ESTIMATE  SIG<=0 */
/*        4      FORGETTING FACTOR   AL<0 OR =>1 */
/*        5      RESPONSE EXCEEDS 1.E10 */

#ifndef KR
void spnlms(float *x, long *n, float *d, float *b, long *l, float *mu, float *sig, float *al, float *px, long *error)
#else
void spnlms(x, n, d, b, l, mu, sig, al, px, error)
long *n, *l, *error;
float *x, *d, *b, *mu, *sig, *al, *px;
#endif
{
    /* Local variables */
    long k, ll;
    float e, tmp;

    if (*l < 0)
    {
	*error = 1;
	return;
    }
    if (*mu <= 0.0 || *mu >= 1.0)
    {
	*error = 2;
	return;
    }
    if (*sig <= 0.0)
    {
	*error = 3;
	return;
    }
    if (*al < 0.0 || *al >= 1.0)
    {
	*error = 4;
	return;
    }

    for (k = 0 ; k < *n ; ++k)
    {
	px[0] = x[k];
	x[k] = 0.0;

	for (ll = 0 ; ll <= *l ; ++ll)
	{
	    x[k] += b[ll] * px[ll];
	}

	if (ABS(x[k]) > BIG)
	{
	    *error = 5;
	    return;
	}

	e = d[k] - x[k];
	*sig = *al * (px[0] * px[0]) + (1.0 - *al) * *sig;
	tmp = *mu * 2 / ((1.0 + (float) *l) * *sig);

	for (ll = 0 ; ll <= *l ; ++ll)
	{
	    b[ll] += tmp * e * px[ll];
	}

	for (ll = *l ; ll >= 1 ; --ll)
	{
	    px[ll] = px[ll - 1];
	}
    }

    *error = 0;
    return;
} /* spnlms */
