#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

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

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

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

  //start in lane 1 (middle lane)
  int lane = 1;

  //have a reference velocity close to speed limit
  double ref_vel = 0.0; //mph

  h.onMessage([&ref_vel, &map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
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

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	//vector<double> next_x_vals;
          	//vector<double> next_y_vals;


          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

          	//staying in the lane with smooth transition between points
          		// previous path (list of xy points) size
            	// be careful not to repeat this variable below again
          	int prev_size = previous_path_x.size();

          		// check for collisions
          	if (prev_size > 0){
          		car_s = end_path_s; //car's "future" position at the end of previous path
          	}

          	bool too_close = false;

          		// find ref_v to use
          	for (int i = 0; i < sensor_fusion.size(); i++){
          		//car is in my lane
          		float d = sensor_fusion[i][6];
          		if (d < (2+4*lane+2) && d > (2+4*lane-2)){
          			double vx = sensor_fusion[i][3];
          			double vy = sensor_fusion[i][4];
          			double check_speed = sqrt(vx*vx+vy*vy);
          			double check_car_s = sensor_fusion[i][5];

          			check_car_s += ((double)prev_size*.02*check_speed); //if using previous points can project s value outwards in time

          			//check s values greater than mine and s gap is smaller than 30 meters (arbitrary value)
          			if ((check_car_s > car_s) && ((check_car_s - car_s) < 30)){
          				// do some logic here, e.g. lower ref velocity so we don't crash into the car in front of us,
          				// could also set the flag to try to change lanes
          				//ref_vel = 29.5; // mph
          				too_close = true;
          				//change blindly to the left lane if there is a car ahead.
          				/*if (lane > 0){
          					lane = 0;
          				}
          				*/

          				// ********* BEHAVIOR PLANNER ******

          				// 1. My FSM states will simply be target lanes, i.e. 0, 1, 2

          				// 2. Calculate cost for each FSM state (each lane)
          				vector<double> costs = {0.0, 0.0, 0.0};
          				vector<double> speed_penalty = {0.0, 0.0, 0.0};
          				double collision_penalty;
          				double total_lane_cost;

          				for (int i = 0; i < sensor_fusion.size(); i++){
          					d = sensor_fusion[i][6];
          					vx = sensor_fusion[i][3];
          					vy = sensor_fusion[i][4];
          					check_speed = sqrt(vx*vx+vy*vy);
          					check_car_s = sensor_fusion[i][5];

          					check_car_s += ((double)prev_size*.02*check_speed);

          					// check for cars within 30 meters in front and behind the car
          					if ( ((check_car_s > car_s) && ((check_car_s - car_s) < 30))
          							|| ((check_car_s < car_s) && ((check_car_s - car_s) > -30)) ){
          						// set collision penalty
          						if (d < 4 && d > 0){
          							// do not set collision penalty for the lane our car is in
          							if (!(d < (2+4*lane+2) && d > (2+4*lane-2))){
          								collision_penalty = 999.0;
          								costs[0] += collision_penalty;
          								cout << "collision in lane 0 detected" << endl;
          								cout << "costs[0] = " << costs[0] << endl;
          							}
          							else {
          								collision_penalty = 99.0;
          								costs[0] += collision_penalty;
          								cout << "collision in lane 0 detected" << endl;
          								cout << "costs[0] = " << costs[0] << endl;
          							}
          						}
          						if (d < 8 && d > 4){
          						    if (!(d < (2+4*lane+2) && d > (2+4*lane-2))){
          						    	collision_penalty = 999.0;
          						    	costs[1] += collision_penalty;
          						    	cout << "collision in lane 1 detected" << endl;
          						    	cout << "costs[1] = " << costs[1] << endl;
          						    }
          						    else {
          								collision_penalty = 99.0;
          								costs[1] += collision_penalty;
          								cout << "collision in lane 1 detected" << endl;
          								cout << "costs[1] = " << costs[1] << endl;
          							}
          						}
          						if (d < 12 && d > 8){
          						    if (!(d < (2+4*lane+2) && d > (2+4*lane-2))){
          						    	collision_penalty = 999.0;
          						    	costs[2] += collision_penalty;
          						    	cout << "collision in lane 2 detected" << endl;
          						    	cout << "costs[2] = " << costs[2] << endl;
          						    }
          						    else {
          								collision_penalty = 99.0;
          								costs[2] += collision_penalty;
          								cout << "collision in lane 2 detected" << endl;
          								cout << "costs[2] = " << costs[2] << endl;
          							}
          						}
          					}
          				}

          				/*
          				for (int i = 0; i < 3; i++){
          					speed_penalty = 0.0; //to be implemented
          					collision_penalty = 0.0; //to be implemented
          					total_lane_cost = speed_penalty + collision_penalty;
          					costs.push_back(total_lane_cost);
          				}
          				*/

          				//consider only neighboring lanes
          				/*
          				if (lane == 0){
          					costs.pop_back();
          				}
          				if (lane == 2){
          					costs.erase(costs.begin());
          				}
          				*/

          				for (int i = 0; i < costs.size(); i++) {
          					cout << "costs[" << i << "] = " << costs[i] << endl;
          				}

          				//choose lane with lowest cost - only from neighboring lanes
          				if (lane == 0){
          					vector<double>::iterator best_cost = min_element(begin(costs), end(costs)-1);
          					int best_idx = distance(begin(costs), best_cost);
          					lane = best_idx;
          				}
          				if (lane == 1){
          				    vector<double>::iterator best_cost = min_element(begin(costs), end(costs));
          				    int best_idx = distance(begin(costs), best_cost);
          				    lane = best_idx;
          				}
          				if (lane == 2){
          				    vector<double>::iterator best_cost = min_element(begin(costs)+1, end(costs));
          				    int best_idx = distance(begin(costs), best_cost);
          				    lane = best_idx;
          				}
          				cout << "selected lane #:" << lane << endl;
          				cout << endl;

          				// ********* EO BEH. PLANNER ******
/*
 * The above needs to be improved in the following way:
 * We need to look at other cars in the lane we want to change into first
 * using similar logic to to seeing if there is a car in front of us or not
 * Going back to Frenet, we can check if there is a car in that lane
 * and then we can check if it's in some gap range of s. If it is,
 * then it's not safe to do that lane change. And then maybe if it's not safe to go left,
 * we can try to go right instead.
 * You would have three states to choose from: keep lane, change lane left, change lane right
 * and e.g. you would never change lane left if you saw there was another car within like 30,
 * or maybe 100 meters in front of you or 50 or so meters behind you.
 * And if you're already in that leftmost lane, you don't want to go off the road.
 * So we need a logic when to shift into each of those states in the FSM.
 * We need a better FSM and a good cost function.
 * Try to look into the future and see what's the best lane to be in 5 seconds or so.
 *
 * Note to self: I think instead of KL, LCL, LCR, it's better to have states which represent
 * target lane - since we have just three target lanes. The I won't have to bother about leaving the road, etc.
 *
 * Recommended flow of project work:
 * 1. Get the car moving
 * 2. Get the car to keep its lane
 * 3. Smooth its path (i.e. no jerk when changing lanes (e.g. with spline)
 * 4. Then you enter the world which looks like the quiz at the end of behavioral planning lesson
 * 	  where probably the best way to deal with it is to have a cost function to decide between maneuvers/lanes
 * 	  The cost function should take into account the cost of being in each lane and choose lowest
 * 5. Predict where the cars will be in the future and what your cost is going to be for being in
 * 	  different states and different lanes in the future. You could use a NB Classifier for that
 * 	  - (but there is no off-the-shelf historical data to train the classifier...)
 *
 * 	  Flow would be:
 * 	  - get the data
 * 	  - try to predict where cars will be in the future
 * 	  - make your own decisions in behavior planning
 * 	  - and then ultimately you generate a trajectory
 * 	  For this project you might work backwards:
 * 	  - build a trajectory and assume no other vehicles are around
 * 	  - start to assume other vehicles are there and build a behavior planner, but don't worry about the future
 * 	  - finally, start predicting where other vehicles will be in the future
 *
 * 	  Cost function design:
 * 	  - Always calculate cost for all three lanes
 * 	  - Output the lowest one, i.e. target lane
 * 	  - Two key cost components: target speed + lateral collision avoidance (those two should do the trick)
 * 	  - Inputs: Vehicle positions & speeds
 * 	    (current, i.e. gap for a lane change will need to be high to avoid collision)
 * 	    (when prediction is added, then the gap can be narrower)
 *
 * 	  Prediction:
 * 	  - Simple prediction assuming the other car is going to continue on its current path
 * 	  - Option a: just extrapolate along current lane
 * 	  - Option b: extrapolate along car's heading (might be tricky in curves and depend on time horizon)
 * 	  - Option c: collect data from the simulator and train a NB Classifier (seems like a lot of effort)
 * 	  - I'll start with option a and see if it's enough.
 * 	  - Food for thought: if another car started changing his lane, we should be able to see that he is not in the middle of the lane
 * 	  - and then we can assume that he is probably going into another lane.
 * 	  - Maybe a function which assigns either one or two potential lanes to the other car would do the trick
 *
 * 	  */

          			}
          		}
          	}

          	/*
          	if(too_close){
          		ref_vel -= .224; //.224 equals roughly to acceleration of 5m/s2
          	}
          	else if (ref_vel < 49.5){
          		ref_vel += .224;
          	}
          	*/

          		// note: the above code will adjust the velocity every cycle
          		// we could be more efficient and adjust the velocity between every path point
          		// down in the path planner

          		// create a list of widely spaced (x,y) waypoints, evenly spaced at 30m
          		// later we will interpolate these waypoints with a spline and fill it in with more points
            vector<double> ptsx;
            vector<double> ptsy;

            	// reference x, y, yaw states
            	// either we will reference the starting point as where the car is or at the previous path's end point
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);

            	// if previous path is almost empty, use the car as a starting reference
            if(prev_size < 2){
            	// use two points that make the path tangent to the car
            	double prev_car_x = car_x - cos(car_yaw);
            	double prev_car_y = car_y - sin(car_yaw);

            	ptsx.push_back(prev_car_x);
            	ptsx.push_back(car_x);

            	ptsy.push_back(prev_car_y);
            	ptsy.push_back(car_y);
            }
            	// use the previous path's end point as starting reference
            else {
            	// redefine reference state as previous path end point
            	ref_x = previous_path_x[prev_size-1];
            	ref_y = previous_path_y[prev_size-1];

            	double ref_x_prev = previous_path_x[prev_size-2];
            	double ref_y_prev = previous_path_y[prev_size-2];
            	ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            	// use two points that make the path tangent to the previous path's end point
            	ptsx.push_back(ref_x_prev);
            	ptsx.push_back(ref_x);

            	ptsy.push_back(ref_y_prev);
            	ptsy.push_back(ref_y);
            }

            // in Frenet, add evenly 30m spaced points ahead of the starting reference
            vector<double> next_wp0 = getXY(car_s+30, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp1 = getXY(car_s+60, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp2 = getXY(car_s+90, (2+4*lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            ptsx.push_back(next_wp0[0]);
            ptsx.push_back(next_wp1[0]);
            ptsx.push_back(next_wp2[0]);

            ptsy.push_back(next_wp0[1]);
            ptsy.push_back(next_wp1[1]);
            ptsy.push_back(next_wp2[1]);

            for (int i = 0; i < ptsx.size(); i++) {
            	// shift car's reference angle to zero degrees
            	double shift_x = ptsx[i] - ref_x;
            	double shift_y = ptsy[i] - ref_y;

            	ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
            	ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
            }

            // create a spline
            tk::spline s;

            // set (x,y) points to the spline
            s.set_points(ptsx, ptsy);

            // define the actual (x,y) points we will use for the planner
            vector<double> next_x_vals;
            vector<double> next_y_vals;

            // start with all the previous path points from the last time
            for (int i = 0; i < previous_path_x.size(); i++){
            	next_x_vals.push_back(previous_path_x[i]);
            	next_y_vals.push_back(previous_path_y[i]);
            }

            // calculate how to break up spline points so that we travel at our desired reference velocity
            double target_x = 30.0; // horizon of 30 meters, at target speed 25 m/s a car would need
            						// more than 1s to reach this point. 0.02s * 50 points = 1s,
            						// i.e. if we generate 50 points they will only reach up to
            						// 25 meters ahead (at target speed) or less (at speed < target speed)
            double target_y = s(target_x);
            double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

            double x_add_on = 0; //starting point in car's coordinate system

            // fill up the rest of our path planner after filling it with previous points,
            // here we will always output 50 points

            for (int i = 1; i <= 50 - previous_path_x.size(); i++){

            	// we should be adding / subtracting ref_vel in this loop
            	// let's try it
            	if(too_close){
            		ref_vel -= .224; //.224 equals roughly to acceleration of 5m/s2
            	}
            	else if (ref_vel < 49.5){
            		ref_vel += .224;
            	}

            	double N = (target_dist/(.02*ref_vel/2.24)); // divided by 2.24 to convert from mph to kmh
            	double x_point = x_add_on + (target_x)/N;
            	double y_point = s(x_point);

            	x_add_on = x_point;

            	double x_ref = x_point;
            	double y_ref = y_point;

            	// rotate back to map coordinates after rotating to car's coordinates earlier
            	x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
            	y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

            	// set the points in relation to reference position
            	x_point += ref_x;
            	y_point += ref_y;

            	next_x_vals.push_back(x_point);
            	next_y_vals.push_back(y_point);
            }

          	//staying in the lane:
            /*
          	double dist_inc = 0.3;
          	for (int i = 0; i < 50; i++){
          		double next_s = car_s + (i+1)*dist_inc;
          		double next_d = 6; //staying in the middle lane (6 meters from the waypoint / middle of the road)
          		vector<double> xy = getXY(next_s, next_d, map_waypoints_s, map_waypoints_x, map_waypoints_y);

          		next_x_vals.push_back(xy[0]);
          		next_y_vals.push_back(xy[1]);
          	}
          	*/

          	// driving in circle
          	/*
            double pos_x;
            double pos_y;
            double angle;
            int path_size = previous_path_x.size();

            for(int i = 0; i < path_size; i++)
            {
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }

            if(path_size == 0)
            {
                pos_x = car_x;
                pos_y = car_y;
                angle = deg2rad(car_yaw);
            }
            else
            {
                pos_x = previous_path_x[path_size-1];
                pos_y = previous_path_y[path_size-1];

                double pos_x2 = previous_path_x[path_size-2];
                double pos_y2 = previous_path_y[path_size-2];
                angle = atan2(pos_y-pos_y2,pos_x-pos_x2);
            }

            double dist_inc = 0.5;
            for(int i = 0; i < 50-path_size; i++)
            {
                next_x_vals.push_back(pos_x+(dist_inc)*cos(angle+(i+1)*(pi()/100)));
                next_y_vals.push_back(pos_y+(dist_inc)*sin(angle+(i+1)*(pi()/100)));
                pos_x += (dist_inc)*cos(angle+(i+1)*(pi()/100));
                pos_y += (dist_inc)*sin(angle+(i+1)*(pi()/100));
            }
            */

          	// *********** eo my code *************

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
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
