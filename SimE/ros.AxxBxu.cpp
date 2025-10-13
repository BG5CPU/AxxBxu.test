#include <my_control/controller.h>
#include <ros/ros.h>
#include <geometry_msgs/PoseStamped.h>
#include <geometry_msgs/TwistStamped.h>
#include <mavros_msgs/AttitudeTarget.h>
#include <mavros_msgs/SetMode.h>
#include <mavros_msgs/CommandBool.h>
#include <mavros_msgs/State.h>
#include <my_control/trajectory_generator.h>
#include <mavros_msgs/ActuatorControl.h>
#include <std_msgs/Float64MultiArray.h>
#include <vector>
#include <Eigen/Dense>
#include <cmath>
#include <fstream>
#include <sensor_msgs/Imu.h>
#include <mav_msgs/eigen_mav_msgs.h>
#include <visualization_msgs/Marker.h>
#include <random>
geometry_msgs::PoseStamped current_pose_leader, current_pose_follower;
geometry_msgs::TwistStamped current_velocity;
geometry_msgs::TwistStamped current_bodyrate;
mavros_msgs::State current_state_leader, current_state_follower1;
mavros_msgs::ActuatorControl control_msg, control_msg_leader, control_msg_follower1;
sensor_msgs::Imu imu_msg;
const double pi = 3.1415926536;
bool collision_free_tracking = false;
bool minimumsnap_check = false;

void poseCallback_leader(const geometry_msgs::PoseStamped::ConstPtr& msg) {current_pose_leader = *msg;}
void poseCallback_follower(const geometry_msgs::PoseStamped::ConstPtr& msg) {current_pose_follower = *msg;}
void velocityCallback(const geometry_msgs::TwistStamped::ConstPtr& msg) {current_velocity = *msg;}
void imuCallback(const sensor_msgs::Imu::ConstPtr& msg) {imu_msg = *msg;}
void stateCallback_leader(const mavros_msgs::State::ConstPtr& msg) {current_state_leader = *msg;}
void stateCallback_follower1(const mavros_msgs::State::ConstPtr& msg) {current_state_follower1 = *msg;}
//  [min, max] 
double uniRand(double min, double max) {
    static std::mt19937 gen(std::random_device{}());  
    std::uniform_real_distribution<double> dist(min, max);
    return dist(gen);
}

static Eigen::Vector3d g_last_euler = Eigen::Vector3d::Zero();
static Eigen::Quaterniond g_last_q(1,0,0,0);
static bool g_inited = false;

Eigen::Vector3d quatToEulerZYX(const Eigen::Quaterniond& q_raw) {
    static bool inited = false;
    static Eigen::Quaterniond q_last(1,0,0,0);
    static Eigen::Vector3d euler_last = Eigen::Vector3d::Zero();

    Eigen::Quaterniond q = q_raw.normalized();
    if (inited && q.dot(q_last) < 0.0) {
        q.coeffs() *= -1.0;
    }
    Eigen::Matrix3d R = q.toRotationMatrix();
    double yaw   = std::atan2(R(1,0), R(0,0));
    double pitch = std::asin(-R(2,0));
    double roll  = std::atan2(R(2,1), R(2,2));
    Eigen::Vector3d euler_new(yaw, pitch, roll);
    if (inited) {
        for (int i=0; i<3; i++) {
            double d = euler_new[i] - euler_last[i];
            if (d >  M_PI) euler_new[i] -= 2*M_PI;
            if (d < -M_PI) euler_new[i] += 2*M_PI;
        }
    }

    q_last = q;
    euler_last = euler_new;
    inited = true;
    return euler_new;
}

class TrajectoryVisualizer {
public:
    TrajectoryVisualizer(){
        // ROS NodeHandle
        ros::NodeHandle nh;
        // Initialize publisher for Marker
        marker_pub_ = nh.advertise<visualization_msgs::Marker>("trajectory_markers", 10);
    }

private:
    ros::Publisher marker_pub_;
};

int main(int argc, char **argv)
{
  setlocale(LC_ALL,"");
  ros::init(argc, argv, "nmpc_zero");
  ros::NodeHandle nh;                  
  ros::NodeHandle pnh("~");            
  ros::Publisher thrust_pub_leader = nh.advertise<mavros_msgs::ActuatorControl>("/HR_leader/mavros/actuator_control",10);
  ros::Publisher LLD_input_pub = nh.advertise<geometry_msgs::PoseStamped>("LLD_input/input_x_y_z", 10);
  ros::Publisher LLD_check_pub = nh.advertise<geometry_msgs::PoseStamped>("LLD_input/attitude_r_p_y", 10);
  ros::Subscriber pose_sub_leader = nh.subscribe("/HR_leader/mavros/local_position/pose", 10, poseCallback_leader);
  ros::Subscriber velocity_sub = nh.subscribe("/HR_leader/mavros/local_position/velocity_local", 10, velocityCallback);
  ros::Subscriber imu_sub = nh.subscribe("/HR_leader/mavros/imu/data",10, imuCallback);
  // offboard check
  ros::Subscriber state_sub_leader = nh.subscribe("/HR_leader/mavros/state", 10, stateCallback_leader);
  ros::spinOnce(); 

  ros::Rate rate(500);

  // std::random_device rd;
  // std::mt19937 gen(rd());  
  // std::uniform_real_distribution<> dis(-0.4, 0.4);  

  while(ros::ok()){
    quadrotor_common::QuadStateEstimate sta_est;
    sta_est.timestamp = ros::Time::now();
    sta_est.coordinate_frame = quadrotor_common::QuadStateEstimate::CoordinateFrame::LOCAL;
    sta_est.position.x() = current_pose_leader.pose.position.x;
    sta_est.position.y() = current_pose_leader.pose.position.y;
    sta_est.position.z() = current_pose_leader.pose.position.z;
    sta_est.velocity.x() = current_velocity.twist.linear.x;
    sta_est.velocity.y() = current_velocity.twist.linear.y;
    sta_est.velocity.z() = current_velocity.twist.linear.z;
    sta_est.orientation.w() = current_pose_leader.pose.orientation.w;
    sta_est.orientation.x() = current_pose_leader.pose.orientation.x;
    sta_est.orientation.y() = current_pose_leader.pose.orientation.y;
    sta_est.orientation.z() = current_pose_leader.pose.orientation.z;
    sta_est.bodyrates.x() = imu_msg.angular_velocity.x;
    sta_est.bodyrates.y() = imu_msg.angular_velocity.y;
    sta_est.bodyrates.z() = imu_msg.angular_velocity.z;
    
    if (true) {
      Eigen::Quaterniond q(
          sta_est.orientation.w(),
          sta_est.orientation.x(),
          sta_est.orientation.y(),
          sta_est.orientation.z()
      );
      Eigen::Vector3d ypr = quatToEulerZYX(q);
      double yaw   = ypr[0];
      double pitch = ypr[1];
      double roll  = ypr[2];
      geometry_msgs::PoseStamped LLD_attitude;
      LLD_attitude.header.stamp = ros::Time::now();
      LLD_attitude.header.frame_id = "map";
      LLD_attitude.pose.position.x = roll;
      LLD_attitude.pose.position.y = pitch;
      LLD_attitude.pose.position.z = yaw;
      LLD_check_pub.publish(LLD_attitude);
      
      //state [Roll, Pitch, Yaw, omega.x, omega.y, omega.z]^T
      Eigen::VectorXd state(6);
      state << roll,
              pitch,
              yaw,
              sta_est.bodyrates.x(),
              sta_est.bodyrates.y(),
              sta_est.bodyrates.z();
      
      //example vk
      Eigen::Matrix<double, 3, 6> vK;
      vK <<-0.8256,    0.4313,   -0.0665,   -0.1294,   -0.0181,   -0.3742,
      -0.2642,    1.4096,   -0.0626,    0.0144,    0.4716,   -0.2499,
      -0.1633,   -0.3253,    2.4733,    0.0767,    0.1194,    1.5078;

      // input = vK * state
      Eigen::Vector3d input = vK * state;
      // data driven controller
      control_msg_leader.controls[0] = input(0); // input_x
      control_msg_leader.controls[1] = input(1); // input_y
      control_msg_leader.controls[2] = input(2); // input_z
      control_msg_leader.controls[3] = 0.8;      

      // record data uniRand(19, 23) 0.1 0.25
      ros::Time t = ros::Time::now();
      double time_now = t.toSec();
      double eu = 0.25; 
      //collect data
      double tau_x = eu * sin(2 * M_PI * 19.8 * time_now + M_PI/7.7);
      double tau_y = eu * sin(2 * M_PI * 20.3 * time_now + M_PI/9.3);
      double tau_z = eu * sin(2 * M_PI * 9.7  * time_now + M_PI/3.2);
      // control_msg_leader.controls[0] = tau_x; // input_x
      // control_msg_leader.controls[1] = tau_y; // input_y
      // control_msg_leader.controls[2] = tau_z; // input_z
      // control_msg_leader.controls[3] = 1.0;

      geometry_msgs::PoseStamped LLD_input;
      LLD_input.header.stamp = ros::Time::now();
      LLD_input.header.frame_id = "map";
      // LLD_input.pose.position.x = tau_x;
      // LLD_input.pose.position.y = tau_y;
      // LLD_input.pose.position.z = tau_z;
      LLD_input.pose.position.x = input(0);
      LLD_input.pose.position.y = input(1);
      LLD_input.pose.position.z = input(2);
      LLD_input_pub.publish(LLD_input);

      ROS_INFO_STREAM("AxxBxu.cpp)"<<std::endl
      << ", Roll: " << roll * 57.63
      << ", Pitch: " << pitch * 57.63
      << ", Yaw: " << yaw * 57.63<< std::endl
      << "tau x: " << control_msg_leader.controls[0] 
      << ", y: " << control_msg_leader.controls[1] 
      << ", z: " << control_msg_leader.controls[2] 
      << ", Thrust: "<< control_msg_leader.controls[3]<<std::endl
      );
      thrust_pub_leader.publish(control_msg_leader);
    } else
    {
      std::cerr << "Error: check flag." << std::endl;
    }
  rate.sleep();
  ros::spinOnce();
  }
  return 0;
}