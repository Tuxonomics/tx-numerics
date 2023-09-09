
// prototype for f: T f( T x[], size_t n )
template <typename T, typename F>
void fd_grad( T *fval, T g[], T x[], T eps[], size_t n, F f )
{
    if ( fval ) {
        *fval = f( x, n );
    }

    T f0 = f( x, n );

    for ( size_t i = 0; i < n; ++i ) {
        x[i] += eps[i];
        g[i] = ( f( x, n ) - f0 ) / eps[i];
        x[i] -= eps[i];
    }
}

template <typename T, typename F>
void fd_grad( double *fval, double g[], double x[], double eps[], size_t n, F f )
{
    fd_grad( fval, g, x, eps, n, f );
}

template <typename T, typename F>
void fd_grad( float *fval, float g[], float x[], float eps[], size_t n, F f )
{
    fd_grad( fval, g, x, eps, n, f );
}



// prototype for f: T f( T x[], size_t n )
template <typename T, typename F>
void cd_grad( T *fval, T g[], T x[], T eps[], size_t n, F f )
{
    if ( fval ) {
        *fval = f( x, n );
    }

    T f0 = T(NAN), f1 = T(NAN), e = T(NAN), eh = T(NAN);

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

template <typename T, typename F>
void cd_grad( double *fval, double g[], double x[], double eps[], size_t n, F f )
{
    cd_grad( fval, g, x, eps, n, f );
}

template <typename T, typename F>
void cd_grad( float *fval, float g[], float x[], float eps[], size_t n, F f )
{
    cd_grad( fval, g, x, eps, n, f );
}

