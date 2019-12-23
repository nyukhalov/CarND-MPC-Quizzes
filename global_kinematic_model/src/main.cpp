// In this quiz you'll implement the global kinematic model.
#include <math.h>
#include <iostream>
#include "Eigen-3.3/Eigen/Core"

using Eigen::VectorXd;

//
// Helper functions
//
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

const double Lf = 2;

// Return the next state.
VectorXd globalKinematic(const VectorXd &state, 
                         const VectorXd &actuators, double dt);

int main() {
  // [x, y, psi, v]
  VectorXd state(4);
  // [delta, a]
  VectorXd actuators(2);

  state << 0, 0, deg2rad(45), 1;
  actuators << deg2rad(5), 1;

  // should be [0.212132, 0.212132, 0.798488, 1.3]
  auto next_state = globalKinematic(state, actuators, 0.3);

  std::cout << next_state << std::endl;
}

VectorXd globalKinematic(const VectorXd &state, 
                         const VectorXd &actuators, double dt) {
  // NOTE: state is [x, y, psi, v] and actuators is [delta, a]

  double heading = state[2];
  double vel = state[3];
  double delta = actuators[0];
  double accel = actuators[1];

  // Create a new vector for the next state.
  VectorXd next_state(state.size());

  // x
  next_state[0] = state[0] + vel*cos(heading)*dt;
  next_state[1] = state[1] + vel*sin(heading)*dt;
  next_state[2] = state[2] + delta*(vel/Lf)*dt;
  next_state[3] = state[3] + accel*dt;

  return next_state;
}
