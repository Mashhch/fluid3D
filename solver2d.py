package fluid;

import javafx.util.Pair;

public class Solver {

    public static Pair<float[][], float[][]> SWAP(float[][] x0,float[][] x) {
        return new Pair<>(x,x0);
    }
    public static void add_source ( int N, float[][] x, float[][] s, float dt )
    {
        for (int i=0 ; i<N+2 ; i++) {
            for (int j=0 ; j<N+2 ; j++ )
                x[i][j] += dt*s[i][j];
        }
    }

    public static void set_bnd ( int N, int b, float [][] x )
    {
        for (int  i=1 ; i<=N ; i++ ) {
            x[0][i] = b==1 ? -x[1][i] : x[1][i];
            x[N+1][i] = b==1 ? -x[N][i] : x[N][i];
            x[i][0] = b==2 ? -x[i][1] : x[i][1];
            x[i][N+1] = b==2 ? -x[i][N] : x[i][N];
        }
        x[0][0] = 0.5f*(x[1][0]+x[0][1]);
        x[0][N+1] = 0.5f*(x[1][N+1]+x[0][N]);
        x[N+1][0] = 0.5f*(x[N][0]+x[N+1][1]);
        x[N+1][N+1] = 0.5f*(x[N][N+1]+x[N+1][N]);
    }

    public static void lin_solve ( int N, int b, float[][] x, float[][] x0, float a, float c )
    {
        for (int k=0 ; k<20 ; k++ ) {
            for (int i=1 ; i<=N ; i++ ) {
                for (int j = 1; j <= N; j++)
                    x[i][j] = (x0[i][j] + a * (x[i - 1][j] + x[i + 1][j] + x[i][j - 1] + x[i][j + 1])) / c;
            }
            set_bnd ( N, b, x );
        }
    }

    public static void diffuse ( int N, int b, float[][] x, float[][] x0, float diff, float dt )
    {
        float a=dt*diff*N*N;
        lin_solve ( N, b, x, x0, a, 1+4*a );
    }

    public static void advect ( int N, int b, float[][] d, float[][] d0, float[][] u, float[][] v, float dt )
    {
        int i0, j0, i1, j1;
        float x, y, s0, t0, s1, t1, dt0;

        dt0 = dt*N;
        for (int i=1 ; i<=N ; i++ ) {
            for (int j = 1; j <= N; j++) {
                x = i - dt0 * u[i][j];
                y = j - dt0 * v[i][j];

                if (x < 0.5f) x = 0.5f;
                if (x > N + 0.5f) x = N + 0.5f;
                i0 = (int) x;
                i1 = i0 + 1;
                if (y < 0.5f) y = 0.5f;
                if (y > N + 0.5f) y = N + 0.5f;
                j0 = (int) y;
                j1 = j0 + 1;
                s1 = x - i0;
                s0 = 1 - s1;
                t1 = y - j0;
                t0 = 1 - t1;
                d[i][j] = s0 * (t0 * d0[i0][j0] + t1 * d0[i0][j1]) +
                        s1 * (t0 * d0[i1][j0] + t1 * d0[i1][j1]);
            }
        }
        set_bnd ( N, b, d );
    }

    public static void project ( int N, float[][] u, float[][] v, float[][] p, float[][] div )
    {
        int i, j;

        for ( i=1 ; i<=N ; i++ ) {
            for ( j=1 ; j<=N ; j++ ) {
                div[i][j] = -0.5f*(u[i+1][j]-u[i-1][j]+v[i][j+1]-v[i][j-1])/N;
                p[i][j] = 0;
            }
        }
        set_bnd ( N, 0, div ); set_bnd ( N, 0, p );
        lin_solve ( N, 0, p, div, 1, 4 );

        for ( i=1 ; i<=N ; i++ ) {
            for (j = 1; j <= N; j++) {
                u[i][j] -= 0.5f * N * (p[i + 1][j] - p[i - 1][j]);
                v[i][j] -= 0.5f * N * (p[i][j + 1] - p[i][j - 1]);
            }
        }
        set_bnd ( N, 1, u ); set_bnd ( N, 2, v );
    }

    public static void dens_step ( int N, float[][] x, float[][] x0, float[][] u, float[][] v, float diff, float dt )
    {
        add_source ( N, x, x0, dt );
        Pair<float[][], float[][]> swap_ = SWAP ( x0, x );
        x0 = swap_.getKey();
        x = swap_.getValue();
        diffuse ( N, 0, x, x0, diff, dt );
        swap_ = SWAP ( x0, x );
        x0 = swap_.getKey();
        x = swap_.getValue();
        advect ( N, 0, x, x0, u, v, dt );
    }

    public static void vel_step ( int N, float[][] u, float[][] v, float[][] u0, float[][] v0, float visc, float dt )
    {
        add_source ( N, u, u0, dt );
        add_source ( N, v, v0, dt );
        Pair<float[][], float[][]> swap_ = SWAP ( u0, u );
        u0 = swap_.getKey();
        u = swap_.getValue();
        diffuse ( N, 1, u, u0, visc, dt );
        swap_= SWAP ( v0, v );
        v0 = swap_.getKey();
        v = swap_.getValue();
        diffuse ( N, 2, v, v0, visc, dt );
        project ( N, u, v, u0, v0 );
        swap_ = SWAP ( u0, u );
        u0 = swap_.getKey();
        u = swap_.getValue();
        swap_ = SWAP ( v0, v );
        v0 = swap_.getKey();
        v = swap_.getValue();
        advect ( N, 1, u, u0, u0, v0, dt );
        advect ( N, 2, v, v0, u0, v0, dt );
        project ( N, u, v, u0, v0 );
    }

}