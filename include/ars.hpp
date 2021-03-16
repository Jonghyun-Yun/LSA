#include "logdensity.h"
#include "lambertW.h"
#include <cmath>
#include <iostream>

using namespace std;

double logsumexp(double u, double v) {
  double rv  = max(u, v) + log(exp(u - max(u, v)) + exp(v - max(u, v)));
  return rv;
}
double logdiffexp(double u, double v) {
  double rv = max(u, v) + log(exp(u - max(u, v)) - exp(v - max(u, v)));
  return rv;
}


struct HullSegment {
public:
  double x;
  double hx;
  double hpx;
  double z_left;
  double z_right;
  double hu_left;
  double hu_right;
  double scum_left;
  double scum_right;
};

class UpperHull {
public:
  int n_segments;
  vector<HullSegment> segments;
  double log_cu; // normalizing constant Z s.t \int e(x) / Z dx = 1 where e(x) = exp(upper.hull(x))
  double xlb;
  double xrb;

  UpperHull() {}; // Default constructor

  UpperHull(vector<double> x, vector<double> hx, vector<double> hpx, double _xlb, double _xrb) {
    // x assumed sorted!

    // No of segments
    n_segments = x.size();

    // Set bounds
    xlb = _xlb;
    xrb = _xrb;

    // Ensure left (right) most point is left (right) of mode
    if (hpx.at(0) < 0) {
      printf("Smallest starting point is right of mode.");
      return;
    }
    if (hpx.at(hpx.size()-1) > 0) {
      printf("Largest starting point is left of mode.");
      return;
    }

    // Initialize segments
    for (int i = 0; i < n_segments; ++i) {

      HullSegment segment;
      segment.x = x.at(i);
      segment.hx = hx.at(i);
      segment.hpx = hpx.at(i);

      segments.push_back(segment);
    }

    // Update z and scum
    update_z();
    update_scum();

  };

  // assign hull values
  void update_z() {

    for (int i = 0; i < n_segments; ++i) {

      HullSegment segment = segments.at(i); // We make a copy and reassign to segments below; inefficient, but safer

      // Assign left z
      if (i == 0) {
        segment.z_left = xlb;
      } else {
        segment.z_left = segments.at(i-1).z_right;
      }

      // Assign right z
      if (i < n_segments - 1) {
        HullSegment next_segment;
        next_segment = segments.at(i+1);
        double num = segment.hx - next_segment.hx + next_segment.hpx*(next_segment.x - segment.x);
        double denom = next_segment.hpx - segment.hpx;
        segment.z_right = segment.x + (num/denom);
      } else {
        segment.z_right = xrb;
      }

      // Assign hull values
      segment.hu_left = segment.hx - segment.hpx*(segment.x - segment.z_left);
      segment.hu_right = segment.hx + segment.hpx*(segment.z_right - segment.x);

      // Reassign segment copy
      segments.at(i) = segment;
    }

  };

  void update_scum() {

    // Compute normalizer
    vector<double> log_cu_vec(n_segments);
    log_cu = 0;
    for (int i = 0; i < n_segments; ++i) {
      HullSegment segment;
      segment = segments.at(i);
      // double cu_i = (1/segment.hpx)*(exp(segment.hu_right) - exp(segment.hu_left));
      double log_cu_i;
      if (segment.hpx > 0) {
        log_cu_i = -log(segment.hpx) + logdiffexp(segment.hu_right, segment.hu_left);
      } else {
        log_cu_i = logdiffexp(segment.hu_left - log(-segment.hpx), segment.hu_right - log(-segment.hpx));
      }

      if (i == 0) {
        log_cu = log_cu_i;
      } else {
        log_cu = logsumexp(log_cu, log_cu_i);
      }
      log_cu_vec.at(i) = log_cu_i;
    }

    // Compute and assign scum
    vector<double> scum_vec(n_segments);
    for (int i = 0; i < n_segments; ++i) {
      if (i == 0) {
        scum_vec.at(0) = exp(log_cu_vec.at(0) - log_cu);
        segments.at(0).scum_left = 0;
        segments.at(0).scum_right = scum_vec.at(0);
      } else {
        scum_vec.at(i) = scum_vec.at(i-1) + exp(log_cu_vec.at(i) - log_cu);
        segments.at(i).scum_left = segments.at(i-1).scum_right;
        segments.at(i).scum_right = scum_vec.at(i);
      }
    }
  };

  // cumulative normalizing constant for each segment
  // scrap from update_scum()
  vector<double> get_log_cu_vec() {

    // Compute normalizer
    vector<double> log_cu_vec(n_segments);
    log_cu = 0;
    for (int i = 0; i < n_segments; ++i) {
      HullSegment segment;
      segment = segments.at(i);
      // double cu_i = (1/segment.hpx)*(exp(segment.hu_right) -
      // exp(segment.hu_left));
      double log_cu_i;
      if (segment.hpx > 0) {
        log_cu_i =
            -log(segment.hpx) + logdiffexp(segment.hu_right, segment.hu_left);
      } else {
        log_cu_i = logdiffexp(segment.hu_left - log(-segment.hpx),
                              segment.hu_right - log(-segment.hpx));
      }

      if (i == 0) {
        log_cu = log_cu_i;
      } else {
        log_cu = logsumexp(log_cu, log_cu_i);
      }
      log_cu_vec.at(i) = log_cu_i;
    }
    return log_cu_vec;
  };

  double sample() {

    double x_sample;
    // double u = R::unif_rand();
    double u = ((double)rand()/(double)RAND_MAX);

    // cout << u << endl;
    double min_scum = segments.at(0).scum_right;
    if (u < min_scum) {
      // Sample is on the left-most segment
      HullSegment s = segments.at(0);

      // Draw a sample from the upper hull

      // x_sample = s.z_right + (1/s.hpx)*log(1 + (s.hpx*exp(log_cu)*(u-s.scum_right))/exp(s.hu_right));
      // ^ This often leads to over/underflow

      if (s.hpx*(u-s.scum_right) > 0) { // LH term might be negative, so we can't compute the log naively
        x_sample = s.z_right + (1/s.hpx)*logsumexp(0, log(s.hpx*(u-s.scum_right)) + log_cu - s.hu_right);
      } else {
        x_sample = s.z_right + (1/s.hpx)*logdiffexp(0, log(-1*s.hpx*(u - s.scum_right)) + log_cu - s.hu_right);
      }

      // cout << "Draw a sample from the upper hull" << endl;

    } else {
      // Determine which segment the sample is on
      int segment_id;
      for (int i = 1; i < n_segments; ++i) {
        if (u > segments.at(i).scum_left && u < segments.at(i).scum_right) {
          segment_id = i;
        }
      }
      HullSegment s = segments.at(segment_id);

      // cout << "Determine which segment the sample is on" << endl;

      // Draw a sample from the upper hull

      // x_sample = s.z_left + (1/s.hpx)*log(1 + (s.hpx*exp(log_cu)*(u-s.scum_left))/exp(s.hu_left));
      // ^ This is prone to over/underflow

      if (s.hpx > 0) { // s.hpx might be negative, so we can't compute the log naively!
        x_sample = s.z_left + (1/s.hpx)*logsumexp(0, log(s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
      } else {
        x_sample = s.z_left + (1/s.hpx)*logdiffexp(0, log(-s.hpx) + log_cu + log(u - s.scum_left) - s.hu_left);
      }
    }

    return x_sample;
  };

  void add_segment(double x, double hx, double hpx) {

    HullSegment segment;
    segment.x = x;
    segment.hx = hx;
    segment.hpx = hpx;

    // Determine segment position
    int iter = 0;
    while(iter < n_segments) {
      if (x < segments.at(iter).x) {
        break;
      }
      iter++;
    }

    // Insert segment
    segments.insert(segments.begin()+iter, segment);
    n_segments = segments.size();

    // Update z and scum
    update_z();
    update_scum();

  };

  // upper hull (log) values
  double get_hu(double x) {

    // Determine which segment x lies on
    int segment_id;
    for (int i = 0; i < n_segments; ++i) {
      if (x > segments.at(i).z_left && x <= segments.at(i).z_right) {
        segment_id = i;
        break;
      }
    }
    HullSegment s = segments.at(segment_id);

    // Get hu
    double hu;
    if (segment_id == 0) {
      hu = s.hu_right - (s.z_right - x)*s.hpx;
    } else {
      hu = s.hu_left + (x - s.z_left)*s.hpx;
    }

    return hu;
  };
};

class LowerHull {
public:
  vector<double> x;
  vector<double> hx;
  vector<double> hpx;
  int n_points;

  LowerHull() {}; // Default constructor

  LowerHull(vector<double> _x, vector<double> _hx, vector<double> _hpx) {
    // x assumed sorted!

    x = _x;
    hx = _hx;
    hpx = _hpx;
    n_points = x.size();

  };

  void add_segment(double new_x, double new_hx, double new_hpx) {

    // Determine point position
    int iter = 0;
    while (iter < n_points) {
      if (new_x < x.at(iter)) {
        break;
      }
      iter++;
    }

    // Assign
    x.insert(x.begin()+iter, new_x);
    hx.insert(hx.begin()+iter, new_hx);
    hpx.insert(hpx.begin()+iter, new_hpx);

    n_points = x.size();
  }

  // lower hull (log) values
  double get_hl(double _x) {

   // Determine point position
   int iter = 0;
    while (iter < n_points) {
      if (_x < x.at(iter)) {
        break;
      }
      iter++;
    }

    double rv;
    if (iter == 0) {
      rv = -numeric_limits<double>::infinity();
    } else if (iter == n_points) {
      rv = -numeric_limits<double>::infinity();
    } else {
      double x_left = x.at(iter-1);
      double x_right = x.at(iter);
      double hx_left = hx.at(iter-1);
      double hx_right = hx.at(iter);
      double d = (hx_right-hx_left)/(x_right-x_left);
      rv = hx_left + (_x-x_left)*d;
    }

    return rv;
  }

};

class ARS {
  // AR sampling using an object of class LogDensity
public:
  UpperHull upper_hull;
  LowerHull lower_hull;
  double xlb;
  double xrb;
  int max_points;

  const static int MAX_REJECTIONS = 500;

  ARS() {}; // Default constructor

  ARS(LogDensity* const log_density, vector<double> x, double xlb, double xrb, int max_points) :
    xlb(xlb),
    xrb(xrb),
    max_points(max_points) {

    // Check bounds
    if (xlb >= xrb) {
      printf("Upper bound is not largetr than lower bound.");
      return;
    }

    // We need at least two starting points
    if (x.size() < 2) {
      printf("At least two starting points required.");
      return;
    }

    // Order x
    sort(x.begin(), x.end());

    // Check points in bounds
    if (x.at(0) < xlb || x.at(x.size()-1) > xrb) {
      printf("Starting point out of bound.");
      return;
    }

    // Evaluate funcs
    vector<double> hx = log_density->h(x);
    vector<double> hpx = log_density->h_prime(x);

    // Try to make starting points valid (if they're invalid)
    int max_tries = 10;
    if (hpx.at(0) <= 0) {
      // Move left-most point left until valid

      double hpx_left = hpx.at(0);
      double x_left = x.at(0);
      int tries = 0;
      while (hpx_left <= 0 && tries < max_tries) {

        if (isfinite(xlb)) {
          x_left -= (x_left-xlb)/2; // Move half-way to the limit
        } else {
          x_left -= pow(2, tries); // Move left by power of 2
        }

        hpx_left = log_density->h_prime(x_left);

        tries++;
      }

      if (tries < max_tries) {
        hpx.at(0) = hpx_left;
        x.at(0) = x_left;
        hx.at(0) = log_density->h(x_left);
      } else {
       printf("Could not find valid lower starting point.");
       return;
      }
    }

    int last_ind = hpx.size() - 1;
    if (hpx.at(last_ind) >= 0) {
      // Move right-most point right until valid

      double hpx_right = hpx.at(last_ind);
      double x_right = x.at(last_ind);
      int tries = 0;
      while (hpx_right >= 0 && tries < max_tries) {

        if (isfinite(xrb)) {
          x_right += (xrb - x_right)/2; // Move half-way to the limit
        } else {
          x_right += pow(2, tries); // Move right by power of 2
        }

        hpx_right = log_density->h_prime(x_right);

        tries++;
      }

      if (tries < max_tries) {
        hpx.at(last_ind) = hpx_right;
        x.at(last_ind) = x_right;
        hx.at(last_ind) = log_density->h(x_right);
      } else {
        // cout << "x candidates: " << x << endl;
        printf("Could not find valid upper starting point.");
        return;
      }
    }

    // Create the hull
    upper_hull = UpperHull(x, hx, hpx, xlb, xrb);
    lower_hull = LowerHull(x, hx, hpx);
  };

  vector<double> sample(unsigned N, LogDensity* const log_density) {

    vector<double> samples;
    int rejections = 0;
    while(samples.size() < N) {

      double x_sample = upper_hull.sample();

      // cout << x_sample << endl;
      // double u = R::runif(0, 1);
      double u = ((double)rand()/(double)RAND_MAX);
      // cout << "hl: " << exp(lower_hull.get_hl(x_sample)) << endl;
      // cout << "hu: " << exp(upper_hull.get_hu(x_sample)) << endl;
      if (u < exp(lower_hull.get_hl(x_sample) - upper_hull.get_hu(x_sample))) {
        // Accept!
        // cout << "Accept!" << endl;
        samples.push_back(x_sample);
        rejections = 0;
      } else {

        // cout << "calculating log_density.. " << endl;

        double hx = log_density->h(x_sample);

        // cout << "hx: " << hx << endl;

        if (u < exp(hx - upper_hull.get_hu(x_sample))) {
          // Accept!
          // cout << "Accept!" << endl;
          samples.push_back(x_sample);
          rejections = 0;
        } else {
          // Reject!
          // cout << "Reject!" << endl;
          rejections++;
        }

        // Add hull segment
        int points = lower_hull.x.size();
        if (points < max_points) {
          double hpx = log_density->h_prime(x_sample);
          upper_hull.add_segment(x_sample, hx, hpx);
          lower_hull.add_segment(x_sample, hx, hpx);
        }
      }

      if (rejections > MAX_REJECTIONS) {
        cout << "Warning: Maximum number of rejections reached. Returning zero sample." << endl;
        samples.push_back(0);
      }

      // checkUserInterrupt();
    }

    return samples;
  };

 vector<double> get_x() {
    return lower_hull.x;
  };

  double get_hu(double x) {
    return upper_hull.get_hu(x);
  };

  double get_hl(double x) {
    return lower_hull.get_hl(x);
  };

  double get_log_cu() {
    return upper_hull.log_cu;
  };

  vector<double> get_log_cu_vec() {
    return upper_hull.get_log_cu_vec();
  };

};

double get_plnp_map(const double &mu, const double &sigma2, const double &y) {
  // Returns the PLN MAP estimate
  double W = lambertW0_CS(sigma2*exp(mu + sigma2*y));
  return -1*W + mu  + sigma2*y;
}

class PLNPosterior {
// Class for sampling from Poisson Log-Normal posterior, i.e.
// p(x_i | y) \propto dpois(y_i | exp(x_i)) * dnorm(x_i | mu_i, sigma_i)
public:

  vector<double> mu;
  Eigen::VectorXd y;
  vector<double> sigma;
  vector<double> segment_x;
  vector<double> log_cu_vec;
  double log_cu;
  vector<ARS> samplers;
  vector<PoissonLogNormal> densities;
  int N;

  PLNPosterior(Eigen::VectorXd &y, vector<double> mu, vector<double> sigma) :
    mu(mu),
    y(y),
    sigma(sigma) {

    // Check arguments
    N = mu.size();
    if (y.size() != N || sigma.size() != N) {
      printf("y, mu, and sigma must be of same length.");
      return;
    }

    // Make the samplers
    for (int i = 0; i < N; ++i) {

      // PLN log density
      PoissonLogNormal log_density(y(i), mu.at(i), sigma.at(i));
      densities.push_back(log_density);

      // MAP estimate
      double x_map = get_plnp_map(mu.at(i), pow(sigma.at(i), 2), y(i));
      if (!isfinite(x_map)) {
        x_map = log(y(i)+1);
      }

      // Starting points
      vector<double> x = {x_map-1, x_map+1};

      // Sampler
      ARS ars(&densities.at(i), x, -numeric_limits<double>::infinity(), numeric_limits<double>::infinity(), 100);
      samplers.push_back(ars);
    }

            // std::cout << "constructed" <<"\n";
  };

  std::vector <std::vector<double> > sample(int M) {

    std::vector <std::vector<double> > out(M, std::vector<double>(N));
    vector<double> smp_i;
    for (int i = 0; i < N; ++i) {
    smp_i = samplers.at(i).sample(M, &densities.at(i));
      for (int j = 0; j < M; ++j) {
        out[j].at(i) = smp_i.at(j);
      }
    }
    return out;
  };
};

vector<double> set_ars_starting(vector<double> segment_x, vector<double> log_cu_vec, double log_cu) {
  double cumsum = 0;
  int n_segments = segment_x.size();
  vector<double> x_init = {0, 0, 0};
  int FIND_switch = 0;
  double subsum;

  for (int i = 0; ((i < n_segments) && (FIND_switch < 3)); ++i) {
    subsum = exp(log_cu_vec.at(i) - log_cu);
    cumsum += subsum;
    std::cout << "x:" << segment_x.at(i) << std::endl;
    std::cout << "cumsum:" << cumsum << std::endl;

    switch(FIND_switch) {
      case 0 :
        if(cumsum >= 0.05) {
          if (i == 0) {
            x_init.at(0) = segment_x.at(i);}
          else {
            x_init.at(0) = segment_x.at(i) - (cumsum - 0.05) / subsum * (segment_x.at(i) - segment_x.at(i-1));
          }
          FIND_switch = 1;
        }
      case 1 :
        if(cumsum >= 0.5) {
          x_init.at(1) = segment_x.at(i) - (cumsum - 0.5) / subsum * (segment_x.at(i) - segment_x.at(i-1));
          FIND_switch = 2;
        }
      case 2 :
        if(cumsum >= 0.95) {
          x_init.at(2) = segment_x.at(i) - (cumsum - 0.95) / subsum * (segment_x.at(i) - segment_x.at(i-1));
          FIND_switch = 99;
          break;
        }
        }
  }
  return x_init;
}
