#pragma once

#include "bm_types.hpp"

void
benchmark_fastscape_sinks(xt::xtensor<bm_index, 1>& stack,
                          xt::xtensor<bm_index, 1>& receivers,
                          xt::xtensor<bm_scalar, 1>& dist2receviers,
                          const xt::xtensor<bm_scalar, 2>& elevation,
                          const xt::xtensor<bool, 2>& active_nodes,
                          bm_scalar dx,
                          bm_scalar dy);

#ifdef ENABLE_RICHDEM
void
benchmark_fastscape_wei2018(xt::xtensor<bm_index, 1>& stack,
                            xt::xtensor<bm_index, 1>& receivers,
                            xt::xtensor<bm_scalar, 1>& dist2receviers,
                            const xt::xtensor<bm_scalar, 2>& elevation,
                            const xt::xtensor<bool, 2>& active_nodes,
                            bm_scalar dx,
                            bm_scalar dy);

void
benchmark_fastscape_zhou2016(xt::xtensor<bm_index, 1>& stack,
                             xt::xtensor<bm_index, 1>& receivers,
                             xt::xtensor<bm_scalar, 1>& dist2receviers,
                             const xt::xtensor<bm_scalar, 2>& elevation,
                             const xt::xtensor<bool, 2>& active_nodes,
                             bm_scalar dx,
                             bm_scalar dy);
#endif
