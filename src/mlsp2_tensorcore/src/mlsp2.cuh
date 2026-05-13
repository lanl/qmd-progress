typedef enum
{
    yes = 0, // do refinement
    no = 1   // do not do refinement
} refine_t;

typedef enum
{
    fp64 = 0,     // uniform double precision
    fp32 = 1,     // uniform single precision
    fp16_fp32 = 2 // accumulate in single, compute in half
} precision_t;


void mlsp2(double *,
           double *,
           double *,
           int,
           int,
           precision_t,
           refine_t);

void mlsp2_alt(
           double *,
           double *,
           double *,
           int,
           int,
           precision_t,
           refine_t);

void mlsp2_alt2(
           double *,
           double *,
           double *,
           int,
           int,
           precision_t,
           refine_t,
           bool);

void mlsp2_full(double *,
                double *,
                int,
                double,
                double,
                precision_t,
                refine_t);

void mlsp2_newton(double *,
                double *,
                int,
                double,
                double,
                double,
                double,
                int,
                precision_t,
                refine_t);
