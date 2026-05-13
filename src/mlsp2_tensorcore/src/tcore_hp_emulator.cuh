void tcoreSPGemmSymm(cublasHandle_t, 
                     const unsigned,
                     const float*, 
                     half*,
                     half*,
                     float*,
                     float*,
                     float*);

void tcoreHPGemmSymm(cublasHandle_t, 
                     const unsigned,
                     const float*, 
                     half*,
                     half*,
                     float*,
                     const float);

void tcoreSPGemmSymmAlt(cublasHandle_t, 
                     const unsigned,
                     const float*, 
                     half*,
                     half*,
                     float*,
                     float*,
                     const float);

// NOTE: first argument overwritten!
void tcoreSPGemmSymmAlt2(cublasHandle_t, 
                     const unsigned,
                     float*, 
                     half*,
                     half*,
                     float*,
                     const float,
                     const float);

void tcoreSPGemmSymm1(cublasHandle_t,
                      const unsigned,
                      const float*,
                      const float*,
                      half*,
                      half*,
                      half*,
                      half*,
                      float*,
                      float*,
                      float*);

void tcoreSPGemmAsymm(cublasHandle_t, 
                     const unsigned,
                     const float*, 
                     half*,
                     half*,
                     const float*, 
                     half*,
                     half*,
                     float*,
                     float*,
                     float*);