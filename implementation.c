#include<stdio.h>
#include<math.h>
#include <time.h>
#define TICK(X) clock_t X = clock()
#define TOCK(X) printf("time %s: %g sec.\n", (#X), (double)(clock() - (X))*1000/CLOCKS_PER_SEC)

//----------------------------------------------------------------------------------//

struct complex_number {   // structure to store and represent complex numbers
    double real;
    double imag;
};

struct complex_number multiply_complex(struct complex_number a, struct complex_number b){  // multiplies 2 complex numbers
    struct complex_number result = {(a.real * b.real - a.imag * b.imag),
                                    (a.real * b.imag + a.imag * b.real)};

    return result;
};

struct complex_number add_complex(struct complex_number a, struct complex_number b){  // adds 2 complex numbers
    struct complex_number result = {(a.real + b.real),
                                    (a.imag + b.imag)};

    return result;
};

//----------------------------------------------------------------------------------//

void dft(int N, struct complex_number x[], struct complex_number X[]) {
    // Calculate Discrete Fourier Transform of signal x[] of length N. Result is X[].

    // DFT[x(n)] = summation(x[n]*exp(-j*2*pi*k*n/N)) ,  n=0 to N-1 for all k 

    for(int k = 0; k < N; k++) {
        struct complex_number summation;
        summation.real = 0;
        summation.imag = 0;

        for(int n = 0; n < N; n++) {
            // calculate every term in summation
            double theta = -2 * M_PI * (k*n) / N;
            struct complex_number Wn = {cos(theta), sin(theta)};

            struct complex_number term = multiply_complex(x[n], Wn);
            summation = add_complex(summation, term);
        }

        X[k] = summation;

    }
}


//----------------------------------------------------------------------------------//

void idft(int N, struct complex_number X[], struct complex_number x[]) {
    // Calculate Inverse DFT of signal X[] of length N. Result is x[].

    for(int n = 0; n < N; n++) {
        struct complex_number summation;
        summation.real = 0;
        summation.imag = 0;

        for(int k = 0; k < N; k++) {
            // calculate every term in summation
            double theta = 2 * M_PI * (k*n) / N;
            struct complex_number Wn = {cos(theta), sin(theta)};

            struct complex_number term = multiply_complex(X[k], Wn);
            summation = add_complex(summation, term);
        }

        summation.real = summation.real / N;
        summation.imag = summation.imag / N;

        x[n] = summation;

    }
}

//----------------------------------------------------------------------------------//

void convolution(int len_x, int len_h, struct complex_number x[], struct complex_number h[], struct complex_number y[]) {
    // Find linear convolution of 2 signals of lengths len_x and len_h. Result will be of length (len_x + len_h - 1).

    int N = len_x + len_h - 1;  // length of result

    struct complex_number x_new[N];
    struct complex_number h_new[N];

    int i;

    for(i = 0; i < N; i++) {  // Zero padding
        if(i < len_x) {
            x_new[i].real = x[i].real;
            x_new[i].imag = x[i].imag;
        }
        else {
            x_new[i].real = 0;
            x_new[i].imag = 0;
        }

        if(i < len_h) {
            h_new[i].real = h[i].real;
            h_new[i].imag = h[i].imag;
        }
        else {
            h_new[i].real = 0;
            h_new[i].imag = 0;
        }
    }

    struct complex_number X[N];
    struct complex_number H[N];
    dft(N, x_new, X);
    dft(N, h_new, H);

    struct complex_number Y[N];

    for (int j= 0; j < N; j++){
        struct complex_number term = multiply_complex(X[j], H[j]);

        Y[j] = term;
    }

    idft(N, Y, y);


    return;
}


//----------------------------------------------------------------------------------//

void fft(struct complex_number x[], int n, struct complex_number tmp[] ) {
  if(n>1) {         /* otherwise, do nothing and return */
    int k,m;    
    struct complex_number z, w, *vo, *ve;
    ve = tmp; vo = tmp+n/2;
    for(k=0; k<n/2; k++) {
      ve[k] = x[2*k];
      vo[k] = x[2*k+1];
    }
    fft( ve, n/2, x );      /* FFT on even-indexed elements of x[] */
    fft( vo, n/2, x );      /* FFT on odd-indexed elements of x[] */
    for(m=0; m<n/2; m++) {
        w.real = cos(2*M_PI*m/(double)n);
        w.imag = -sin(2*M_PI*m/(double)n);
        z.real = w.real*vo[m].real - w.imag*vo[m].imag; /* Re(w*vo[m]) */
        z.imag = w.real*vo[m].imag + w.imag*vo[m].real; /* Im(w*vo[m]) */
        x[m].real = ve[m].real + z.real;
        x[m].imag = ve[m].imag + z.imag;
        x[m+n/2].real = ve[m].real - z.real;
        x[m+n/2].imag = ve[m].imag - z.imag;
    }
  }
  return;
}

//----------------------------------------------------------------------------------//


int main()
{

    int selector, N, function;

    printf("Press 1 for DFT, 2 for IDFT, 3 for Convolution, 4 for FFT and 5 for Excecution time:");
    scanf("%d", &function);

    if (function == 1){

        printf( "Enter the length of signal i.e. N = ");
        scanf("%d", &N);

        struct complex_number x[N];
        struct complex_number X[N];

        int i;

        printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
        scanf("%d", &selector);

        printf( "Enter the values of x[n] : \n");

        if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
            for(i = 0; i < N; i++) {
                scanf("%lf", &x[i].real);
                x[i].imag = 0;
            }
        }
        else{
            for(i = 0; i < N; i++) {
                printf("Real: ");
                scanf("%lf", &x[i].real);

                printf("Imaginary: ");
                scanf("%lf", &x[i].imag);
            }
        }

        //------------------- DFT ---------------------//

        dft(N, x, X);

        printf("DFT result, X[k]: \n");
        for(i = 0; i < N; i++) {
            printf("X[%d] = %lf + j%lf \n", i, X[i].real, X[i].imag);

        }
    }
    else if (function == 2){
        printf( "Enter the length of signal i.e. N = ");
        scanf("%d", &N);

        struct complex_number x[N];
        struct complex_number X[N];
       
        int i;

        printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
        scanf("%d", &selector);
 
        printf( "Enter the values of x[n] : \n");

        if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
            for(i = 0; i < N; i++) {
                scanf("%lf", &X[i].real);
                X[i].imag = 0;
            }
        }
        else{
            for(i = 0; i < N; i++) {
                printf("Real: ");
                scanf("%lf", &X[i].real);

                printf("Imaginary: ");
                scanf("%lf", &X[i].imag);
            }
        }

        //------------------- IDFT --------------------//

        idft(N, X, x);

        printf("IDFT result:, x[n]: \n");
        for(i = 0; i < N; i++) {
            printf("x[%d] = %lf + j%lf \n", i, x[i].real, x[i].imag);

        }
    }else if (function == 3){
        int x_length, h_length;
        printf( "Enter the length of signal x i.e. len_x = ");
        scanf("%d", &x_length);
        printf( "Enter the length of signal h i.e. len_h = ");
        scanf("%d", &h_length);

        struct complex_number x[x_length];
        struct complex_number h[h_length];

        int output_length = x_length + h_length - 1;  // length of output
        struct complex_number y[output_length];

        int i;

        printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
        scanf("%d", &selector);

        if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
            printf( "Enter the values of x[n] : \n");
            for(i = 0; i < x_length; i++) {
                scanf("%lf", &x[i].real);
                x[i].imag = 0;
            }
            printf( "Enter the values of h[n] : \n");
            for (int i = 0; i < h_length; i++){
                scanf("%lf", &h[i].real);
                h[i].imag = 0;
            }
        }
        else{
            printf( "Enter the values of x[n] : \n");
            for(i = 0; i < x_length; i++) {
                printf("Real: ");
                scanf("%lf", &x[i].real);

                printf("Imaginary: ");
                scanf("%lf", &x[i].imag);
            }
            printf( "Enter the values of h[n] : \n");
            for(i = 0; i < h_length; i++) {
                printf("Real: ");
                scanf("%lf", &h[i].real);

                printf("Imaginary: ");
                scanf("%lf", &h[i].imag);
            }
        }

        convolution(x_length, h_length, x, h, y);

        printf("CONVOLUTION result:, y[n]: \n");
        for(i = 0; i < output_length; i++) {
            printf("y[%d] = %lf + j%lf \n", i, y[i].real, y[i].imag);

        }
    }else if (function == 4 ){
        printf( "Enter the length of signal i.e. N = ");
        scanf("%d", &N);

        struct complex_number x[N];
        struct complex_number temp[N];

        int i;

        printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
        scanf("%d", &selector);

        printf( "Enter the values of x[n] : \n");

        if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
            for(i = 0; i < N; i++) {
                scanf("%lf", &x[i].real);
                x[i].imag = 0;
            }
        }
        else{
            for(i = 0; i < N; i++) {
                printf("Real: ");
                scanf("%lf", &x[i].real);

                printf("Imaginary: ");
                scanf("%lf", &x[i].imag);
            }
        }

        //------------------- DFT ---------------------//

        fft(x, N, temp);

        printf("FFT result, X[k]: \n");
        for(i = 0; i < N; i++) {
            printf("X[%d] = %lf + j%lf \n", i, x[i].real, x[i].imag);

        }
    }else if (function == 5){
        printf( "Enter the length of signal i.e. N = ");
        scanf("%d", &N);

        struct complex_number x[N];
        struct complex_number X[N];
        struct complex_number temp[N];

        int i;

        printf( "Is x[n] a real valued signal? (1: Yes, 0: No): ");
        scanf("%d", &selector);

        printf( "Enter the values of x[n] : \n");

        if (selector == 1){  // If yes (real valued), accept only real values for ease of use.
            for(i = 0; i < N; i++) {
                scanf("%lf", &x[i].real);
                x[i].imag = 0;
            }
        }
        else{
            for(i = 0; i < N; i++) {
                printf("Real: ");
                scanf("%lf", &x[i].real);

                printf("Imaginary: ");
                scanf("%lf", &x[i].imag);
            }
        }

        TICK(TIME_FFT);
        dft(N, x, X);
        TOCK(TIME_FFT);

        TICK(TIME_DFT);
        fft(x, N, temp);
        TOCK(TIME_DFT);

    }

    return 0;

}
