#include <float.h>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"

#include "car.hpp"
#include "jmt.hpp"
#include "json.hpp"
#include "spline.h"

using namespace std;
using namespace Eigen;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2) {
  return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x,
                    vector<double> maps_y) {
  double closestLen = 100000;  // large number
  int closestWaypoint = 0;

  for (int i = 0; i < maps_x.size(); i++) {
    double map_x = maps_x[i];
    double map_y = maps_y[i];
    double dist = distance(x, y, map_x, map_y);
    if (dist < closestLen) {
      closestLen = dist;
      closestWaypoint = i;
    }
  }

  return closestWaypoint;
}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x,
                 vector<double> maps_y) {
  int closestWaypoint = ClosestWaypoint(x, y, maps_x, maps_y);

  double map_x = maps_x[closestWaypoint];
  double map_y = maps_y[closestWaypoint];

  double heading = atan2((map_y - y), (map_x - x));

  double angle = abs(theta - heading);

  if (angle > pi() / 4) {
    closestWaypoint++;
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta,
                         vector<double> maps_x, vector<double> maps_y) {
  int next_wp = NextWaypoint(x, y, theta, maps_x, maps_y);

  int prev_wp;
  prev_wp = next_wp - 1;
  if (next_wp == 0) {
    prev_wp = maps_x.size() - 1;
  }

  double n_x = maps_x[next_wp] - maps_x[prev_wp];
  double n_y = maps_y[next_wp] - maps_y[prev_wp];
  double x_x = x - maps_x[prev_wp];
  double x_y = y - maps_y[prev_wp];

  // find the projection of x onto n
  double proj_norm = (x_x * n_x + x_y * n_y) / (n_x * n_x + n_y * n_y);
  double proj_x = proj_norm * n_x;
  double proj_y = proj_norm * n_y;

  double frenet_d = distance(x_x, x_y, proj_x, proj_y);

  // see if d value is positive or negative by comparing it to a center point

  double center_x = 1000 - maps_x[prev_wp];
  double center_y = 2000 - maps_y[prev_wp];
  double centerToPos = distance(center_x, center_y, x_x, x_y);
  double centerToRef = distance(center_x, center_y, proj_x, proj_y);

  if (centerToPos <= centerToRef) {
    frenet_d *= -1;
  }

  // calculate s value
  double frenet_s = 0;
  for (int i = 0; i < prev_wp; i++) {
    frenet_s += distance(maps_x[i], maps_y[i], maps_x[i + 1], maps_y[i + 1]);
  }

  frenet_s += distance(0, 0, proj_x, proj_y);

  return {frenet_s, frenet_d};
}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s,
                     vector<double> maps_x, vector<double> maps_y) {
  int prev_wp = -1;

  while (s > maps_s[prev_wp + 1] && (prev_wp < (int)(maps_s.size() - 1))) {
    prev_wp++;
  }

  int wp2 = (prev_wp + 1) % maps_x.size();

  double heading =
      atan2((maps_y[wp2] - maps_y[prev_wp]), (maps_x[wp2] - maps_x[prev_wp]));
  // the x,y,s along the segment
  double seg_s = (s - maps_s[prev_wp]);

  double seg_x = maps_x[prev_wp] + seg_s * cos(heading);
  double seg_y = maps_y[prev_wp] + seg_s * sin(heading);

  double perp_heading = heading - pi() / 2;

  double x = seg_x + d * cos(perp_heading);
  double y = seg_y + d * sin(perp_heading);

  return {x, y};
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  int lane = 1;

  h.onMessage([&map_waypoints_x, &map_waypoints_y, &map_waypoints_s,
               &map_waypoints_dx, &map_waypoints_dy,
               &lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                      uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    // auto sdata = string(data).substr(0, length);
    // cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {
      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

          // Main car's localization Data
          double car_x = j[1]["x"];
          double car_y = j[1]["y"];
          double car_s = j[1]["s"];
          double car_d = j[1]["d"];
          double car_yaw = j[1]["yaw"];
          double car_speed = j[1]["speed"];

          // Previous path data given to the Planner
          auto previous_path_x = j[1]["previous_path_x"];
          auto previous_path_y = j[1]["previous_path_y"];

          // Previous path's end s and d values
          double end_path_s = j[1]["end_path_s"];
          double end_path_d = j[1]["end_path_d"];

          // Sensor Fusion Data, a list of all other cars on the same side of
          // the road.
          auto sensor_fusion = j[1]["sensor_fusion"];

          for (auto &sensor_data : sensor_fusion) {
            int id = sensor_data[0];
            double x = sensor_data[1];
            double y = sensor_data[2];
            double vx = sensor_data[3];
            double vy = sensor_data[4];
            double s = sensor_data[5];
            double d = sensor_data[6];
          }

          int prev_size = previous_path_x.size();
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);
          const double max_speed = 49. * 1.609344 * 1000 / 3600;
          double target_velc = max_speed;

          vector<double> ptsx;
          vector<double> ptsy;

          // car_speed in m/s
          car_speed *= 1.609344 * 1000 / 3600;

          if (prev_size < 2) {
            double prev_car_x = car_x - cos(ref_yaw);
            double prev_car_y = car_y - sin(ref_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);
            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          } else {
            ref_x = previous_path_x[prev_size - 1];
            double ref_x_prev = previous_path_x[prev_size - 2];
            ref_y = previous_path_y[prev_size - 1];
            double ref_y_prev = previous_path_y[prev_size - 2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);
            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);

            car_s = end_path_s;

            double current_speed =
                (sqrt((ref_x - ref_x_prev) * (ref_x - ref_x_prev) +
                      (ref_y - ref_y_prev) * (ref_y - ref_y_prev)) /
                 .02);

            car_speed = current_speed;
          }

          // find the target velocity from the closest car
          double closest_dist = DBL_MAX;
          bool change_lanes = false;

          for (int i = 0; i < sensor_fusion.size(); i++) {
            float d = sensor_fusion[i][6];
            if (d < (2 + 4 * lane + 2) && d > (2 + 4 * lane - 2)) {
              double vx = sensor_fusion[i][3];
              double vy = sensor_fusion[i][4];
              double check_speed = sqrt(vx * vx + vy * vy);
              double check_car_s = sensor_fusion[i][5];
              // make prediction
              check_car_s += prev_size * .02 * check_speed;

              if ((check_car_s > car_s) && ((check_car_s - car_s) < 30) &&
                  ((check_car_s - car_s) < closest_dist)) {
                closest_dist = check_car_s - car_s;

                change_lanes = true;
                
                if ((check_car_s - car_s) > 20) {
                  target_velc = check_speed;
                } else {
                  target_velc = check_speed * 0.8;
                }
              }
            }
          }

          if (!change_lanes) target_velc = max_speed;

          bool changed_lanes = false;

          if (change_lanes) {
            // try the left lane
            if (lane != 0 && !changed_lanes) {
              bool lane_safe = true;
              for (int i = 0; i < sensor_fusion.size(); i++) {
                float d = sensor_fusion[i][6];
                if (d < (2 + 4 * (lane - 1) + 2) &&
                    d > (2 + 4 * (lane - 1) - 2)) {
                  double vx = sensor_fusion[i][3];
                  double vy = sensor_fusion[i][4];
                  double check_speed = sqrt(vx * vx + vy * vy);

                  double check_car_s = sensor_fusion[i][5];
                  check_car_s += prev_size * .02 * check_speed;
                  double dist_s = check_car_s - car_s;
                  if (dist_s < 20 && dist_s > -20) {
                    lane_safe = false;
                  }
                }
              }
              if (lane_safe) {
                changed_lanes = true;
                lane -= 1;
              }
            }

            // try the right lane
            if (lane != 2 && !changed_lanes) {
              bool lane_safe = true;
              for (int i = 0; i < sensor_fusion.size(); i++) {
                // car is in right lane
                float d = sensor_fusion[i][6];
                if (d < (2 + 4 * (lane + 1) + 2) &&
                    d > (2 + 4 * (lane + 1) - 2)) {
                  double vx = sensor_fusion[i][3];
                  double vy = sensor_fusion[i][4];
                  double check_speed = sqrt(vx * vx + vy * vy);

                  double check_car_s = sensor_fusion[i][5];
                  check_car_s += prev_size * .02 * check_speed;
                  double dist_s = check_car_s - car_s;
                  cout << "dist_s: " << dist_s << endl;
                  if (dist_s < 20 && dist_s > -20) {
                    lane_safe = false;
                  }
                }
              }
              if (lane_safe) {
                changed_lanes = true;
                lane += 1;
              }
            }
          }

          vector<double> next_wp0 =
              getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s,
                    map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 =
              getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s,
                    map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 =
              getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s,
                    map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);
          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          // transform the way points to local coordinates
          Matrix3d m_rot;
          m_rot << cos(ref_yaw), -sin(ref_yaw), 0., sin(ref_yaw), cos(ref_yaw),
              0, 0., 0., 1.;
          Matrix3d m_trans;
          m_trans << 1., 0., ref_x, 0., 1., ref_y, 0., 0., 1.;
          Matrix3d m_transform = m_trans * m_rot;
          Matrix3d m_coord = m_transform.inverse();
          for (int i = 0; i < ptsx.size(); i++) {
            Vector3d v;
            v << ptsx[i], ptsy[i], 1.;
            Vector3d p = m_coord * v;
            ptsx[i] = p(0);
            ptsy[i] = p(1);
          }

          // drive
          json msgJson;
          vector<double> next_x_vals;
          vector<double> next_y_vals;

          tk::spline spline;
          spline.set_points(ptsx, ptsy);

          for (int i = 0; i < previous_path_x.size(); i++) {
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double target_x = 30.0;
          double target_y = spline(target_x);
          double target_dist =
              sqrt((target_x) * (target_x) + (target_y) * (target_y));

          double x_add_on = 0;

          for (int i = 1; i <= 50 - previous_path_x.size(); i++) {
            if (car_speed < target_velc) {
              car_speed += target_velc * 0.02 * 0.2;
            } else {
              car_speed -= target_velc * 0.02 * 0.2;
            }

            double N = (target_dist / (.02 * car_speed));
            double x_point = x_add_on + (target_x) / N;
            double y_point = spline(x_point);

            x_add_on = x_point;

            double x_ref = x_point;
            double y_ref = y_point;

            Vector3d v;
            v << x_point, y_point, 1.;
            Vector3d p = m_transform * v;
            next_x_vals.push_back(p(0));
            next_y_vals.push_back(p(1));
          }

          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;
          auto msg = "42[\"control\"," + msgJson.dump() + "]";
          // this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
