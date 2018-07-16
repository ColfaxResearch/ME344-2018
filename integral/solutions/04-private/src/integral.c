#include <math.h>

float IntegrateMyFunction(int const n, float const a, float  const b) {
  
  // Running sum of the integral
  float I = 0.0f;

  // Integration interval
  float const dx = (b-a)/float(n);

  // Loop through the integration range
#pragma omp parallel 
  {

    // Private to the current thread
    float I_partial = 0.0f;

#pragma omp for
    for (int i = 0; i < n; i++) {

      // Midpoint of the integration interval
      float const x = a + dx*(float(i) + 0.5f);

      // Function value at the midpoint
      float const f = 1.0f/sqrtf(x);
    
      // Incrementing the running partial sum
      I_partial += f;

    }
    
    // Reduction (aggregation) into the shared counter
#pragma omp atomic
    I += I_partial;

  }

  // Scale according to the integration interval
  I *= dx;

  return I;
  
}
