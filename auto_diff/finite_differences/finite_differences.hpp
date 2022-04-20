
// prototype for f: double f( double x[], size_t n )
template <typename F>
void fd_grad( double *fval, double g[], double x[], double eps[], size_t n, F f )
{
    if ( fval ) {
        *fval = f( x, n );
    }

    double f0 = f( x, n );

    for ( size_t i = 0; i < n; ++i ) {
        x[i] += eps[i];
        g[i] = ( f( x, n ) - f0 ) / eps[i];
        x[i] -= eps[i];
    }
}



// prototype for f: double f( double x[], size_t n )
template <typename F>
void cd_grad( double *fval, double g[], double x[], double eps[], size_t n, F f )
{
    if ( fval ) {
        *fval = f( x, n );
    }

    double f0 = NAN, f1 = NAN, e = NAN, eh = NAN;

    for ( size_t i = 0; i < n; ++i ) {
        e  = eps[i];
        eh = e/2.0;

        x[i] += eh;
        f0 = f( x, n );
        x[i] -= e;
        f1 = f( x, n );
        x[i] += eh;

        g[i] = ( f1 - f0 ) / e;
    }
}


