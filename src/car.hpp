#ifndef PATH_PLANNING_CAR
#define PATH_PLANNING_CAR

using namespace std;

class EgoCar {

};

class Car {
 public:
  int id_;
  double s_;
  double d_;
  double v_;

  Car(int id, double s, double d, double v)
      : id_(id), s_(s), d_(d), v_(v) {}
  
  double pred_s(double t = 1) {
    return s_ + v_ * t;
  }

  vector<double> state_in(double t = 1) {
    return vector<double>();
  }
};

#endif