#ifndef PI
#define PI 3.14159265358979323846
#endif

/* #ifndef my_dim */
/* #define my_dim 5 */
/* #endif */

#ifndef REPEAT
#define REPEAT 1
#endif

#ifndef DIST_HPP //  This is defined in "dist.hpp" so everything will be skipped
                 //  until #endif is encountered.
#define DIST_HPP
#include "dist.hpp"
#endif

#ifndef DIST2_HPP
#define DIST2_HPP
#include "dist2.hpp"
#endif

#ifndef GLOBAL_HEADER
#define GLOBAL_HEADER
#endif

/* #ifndef nwarm */
// #define nwarm 2057500 // samc dim 11 18sec
// #define nwarm 1600000 // samc dim 11 11.7 sec
// #define nwarm 644444 // samc dim 3
// #define nwarm 200000
/* #define nwarm 830000 // samc dim 5 */
/* #endif */

/* #ifndef niter */
// #define nwarm 8230000 // samc dim 11 18sec
// #define niter 5200000 // samc dim 11 11.7 sec
// #define niter 2577778 // samc dim 3
// #define niter 800000
/* #define niter 3320000 // samc dim 5 */
/* #endif */

/* #ifndef my_eps */
/* #define my_eps 0.971 */
// #define my_eps 0.9 // dim 3
/* #define my_eps 0.9 // dim 3 */
/* #endif */

/* #ifndef my_L */
// #define my_L 1 // dim 3
/* #define my_L 3 // dim 3 */
/* #endif */
