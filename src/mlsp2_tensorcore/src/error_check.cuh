#include <cufft.h>
#include <cusolverDn.h>
#include <cublas_v2.h>

#define CUDA_CHECK_ERR(val) check((val), #val, __FILE__, __LINE__)
#define CUFFT_CHECK_ERR(ans)                    \
    {                                           \
        cufft_check((ans), __FILE__, __LINE__); \
    }
#define CUSOLVER_CHECK_ERR(ans)                    \
    {                                              \
        cusolver_check((ans), __FILE__, __LINE__); \
    }
#define CUBLAS_CHECK_ERR(ans)                    \
    {                                            \
        cublas_check((ans), __FILE__, __LINE__); \
    }

inline void cublas_check(int code, const char *file, int line, bool abort = true)
{
    if (code != CUBLAS_STATUS_SUCCESS)
    {
        fprintf(stderr, "CUBLAS_CHECK_ERR: %d %s %d\n", code, file, line);
        if (abort)
            exit(code);
    }
}

inline void cusolver_check(int code, const char *file, int line, bool abort = true)
{
    if (code != CUSOLVER_STATUS_SUCCESS)
    {
        fprintf(stderr, "CUSOLVER_CHECK_ERR: %d %s %d\n", code, file, line);
        if (abort)
            exit(code);
    }
}

inline void cufft_check(int code, const char *file, int line, bool abort = true)
{
    if (code != CUFFT_SUCCESS)
    {
        fprintf(stderr, "CUFFT_CHECK_ERR: %d %s %d\n", code, file, line);
        if (abort)
            exit(code);
    }
}

template <typename T>
void check(T err, const char *const func, const char *const file,
           const int line)
{
    if (err != cudaSuccess)
    {
        std::cerr << "CUDA Runtime Error at: " << file << ":" << line
                  << std::endl;
        std::cerr << cudaGetErrorString(err) << " " << func << std::endl;
    }
}
