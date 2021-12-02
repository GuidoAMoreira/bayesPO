#include "progressBar.h"

void MultiProgressBar::initializer() {
  if (starter) {
    if (warmupPercentage > 0) {
      for (int i = 0; i < warmupPercentage * _max_ticks; i++)
        REprintf(" ");
      REprintf("Warmup ends here\n");
      for (int i = 0; i < warmupPercentage * _max_ticks; i++)
        REprintf(" ");
      REprintf("\\/\n");
    }
    REprintf("    Sampling 10   20   30   40   50   60   70   80   90   100%\n");
    for (int c = 0; c < chains; c++)
      REprintf("Chain %i: [----|----|----|----|----|----|----|----|----|----|\n", c + 1);
    flush_console();
  }

  starter = false;
}
