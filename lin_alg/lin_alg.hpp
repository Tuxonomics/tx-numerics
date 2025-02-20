
void mm_naive(
    double *c, double *a, double *b,
    int n, int m, int p
) {
#define A(i,j) a[i + n*j]
#define B(i,j) b[i + m*j]
#define C(i,j) c[i + n*j]

    for ( int i = 0; i < n; ++i )
    {
        for ( int k = 0; k < p; ++k )
        {
            C(i,k) = 0.0;
            for ( int j = 0; j < m; ++j )
            {
               C(i,k) += A(i,j) * B(j,k);
            }
        }
    }

#undef A
#undef B
#undef C
}


void mm_linear(
    double *c, double *a, double *b,
    int n, int m, int p
) {
#define A(i,j) a[i + n*j]
#define B(i,j) b[i + m*j]
#define C(i,j) c[i + n*j]

    int i = 0;

    for ( int k = 0; k < p; ++k )
    {
        for ( int j = 0; j < m; ++j )
        {
            for ( i = 0; i < n; ++i )
            {
                if ( j == 0 )
                {
                    C(i,k) = 0.0;
                }
                C(i,k) += A(i,j) * B(j,k);
            }
        }
    }

#undef A
#undef B
#undef C
}


void mm_linear_x4(
    double *c, double *a, double *b,
    int n, int m, int p
) {
#define A(i,j) a[i + n*j]
#define B(i,j) b[i + m*j]
#define C(i,j) c[i + n*j]

    double *cc = c;
    double *aa = a;
    double *bb = b;

    double *cend = c + n*p;
    double *aend = a + n*m;
    double *bend = b + m*p;

    int n4 = n - (n % 4);

    int i = 0;

    for ( int k = 0; k < p; ++k )
    {
        for ( int j = 0; j < m; ++j )
        {
            for ( i = 0; i < n4; i += 4 )
            {
                if ( j == 0)
                {
                    C(i,  k) = 0.0;
                    C(i+1,k) = 0.0;
                    C(i+2,k) = 0.0;
                    C(i+3,k) = 0.0;
                }

                C(i,k)   += A(i,j)   * B(j,k);
                C(i+1,k) += A(i+1,j) * B(j,k);
                C(i+2,k) += A(i+2,j) * B(j,k);
                C(i+3,k) += A(i+3,j) * B(j,k);
            }
            for ( ; i < n; i += 1 )
            {
                if ( j == 0 )
                {
                    C(i,k) = 0.0;
                }
                C(i,k) += A(i,j) * B(j,k);
            }
        }
    }

#undef A
#undef B
#undef C
}


void mm_linear_x4_ptr(
    double *c, double *a, double *b,
    int n, int m, int p
) {
    double *cc = c;
    double *aa = a;
    double *bb = b;

    int n4 = n - (n % 4);

    int i = 0;

    for ( int k = 0; k < p; ++k )
    {
        // start at A(0,0) for every new k
        aa = a;

        for ( int j = 0; j < m; ++j )
        {
            // k-th column of C
            cc = c + k*n;

            for ( i = 0; i < n4; i += 4 )
            {
                if ( j == 0)
                {
                    cc[0] = 0.0;
                    cc[1] = 0.0;
                    cc[2] = 0.0;
                    cc[3] = 0.0;
                }

                cc[0] += aa[0] * bb[0];
                cc[1] += aa[1] * bb[0];
                cc[2] += aa[2] * bb[0];
                cc[3] += aa[3] * bb[0];

                aa += 4;
                cc += 4;
            }

            for ( ; i < n; i += 1 )
            {
                if ( j == 0 )
                {
                    cc[0] = 0.0;
                }
                cc[0] += aa[0] * bb[0];

                aa += 1;
            }

            bb += 1;
        }

    }
}




