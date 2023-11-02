#ifndef WIREBENDING_ONCE
#define WIREBENDING_ONCE
#pragma once

#include <Eigen/Geometry>
#include <vector>
#include <igl/readPLY.h>
#include <igl/opengl/glfw/Viewer.h> 
#include <igl/ray_mesh_intersect.h>
#include <igl/read_triangle_mesh.h>
#include <igl/copyleft/cgal/intersect_other.h>
#include <cstdlib>
#include <utility>
#include <iostream>
#include <cmath>


class Wire {
public:
    Wire();

    void initialize(); 
    void modifyPoint(int index, double x, double y, double z); 
    Eigen::MatrixXd points;

private:
    
};

class Bending {
public:
    Bending(Eigen::MatrixXd& V1_, Eigen::MatrixXd& V2_, Eigen::MatrixXd& V3_,
        Eigen::MatrixXi& F1_, Eigen::MatrixXi& F2_, Eigen::MatrixXi& F3_,
        Eigen::MatrixXd& wire_, std::vector<double>& distances_, std::vector<double>& angles_, int& anglesnum_, std::vector<double>& angles2_);


    double calculateAngle(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3);
    bool blenderback_place();//回撤完成位置
    bool blenderback_place2();
    bool blenderdirchange();//下一轮位置是否在完成位置范围内？大于小于需要改变方向否 结束后统一方向
    Eigen::RowVector3d compute_rotation_point(const Eigen::RowVector3d& initial, const Eigen::RowVector3d& rotated, const Eigen::RowVector3d& initial2);//线旋转
    void blenderback();//pin 下+回撤
    void wireforward();//线走+pin 上
    void wirebending();//弯曲
    void drawpic();//绘制二维坐标轴
    
    
    //bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);//控制输入
    void drawAnimation();

    Eigen::MatrixXd V1, V2, V3;
    Eigen::MatrixXi F1, F2, F3;
    Eigen::MatrixXd w;
    std::vector<double> distances;
    std::vector<double> angles;
    std::vector<double>& angles2;
    igl::opengl::glfw::Viewer viewer;
    int step = 0;
    int anglesnum;
    Eigen::RowVector3d a;//判断旋转角度是否完成
    bool blenderdirchanged = 0;//是否未在范围内并且需要改变方向才能进入范围
    int updownstep = 0;//上下步骤
    int updownstepcount = 0;//上下多少
    double z_offset = 0.5;//上下速度
    double theta = 0.03;// 控制旋转角度和方向的变化
    Eigen::RowVector3d V2point_center = (V2.row(112) + V2.row(205) + V2.row(173) + V2.row(148)) / 4;
    Eigen::RowVector3d jichupoint;
    Eigen::RowVector3d center3;
    int already_length[100] = { 0 };
    std::vector<igl::Hit> hits;//线动
    Eigen::MatrixXi IF;
    bool intersection = 0;//检测碰撞
    Eigen::RowVector3d pos;//碰撞点
    Eigen::RowVector3d pos2;//碰撞点之前位置
    double intersection_1, intersection_2, intersection_3, intersection_4, intersection_5;  //ef 需要 double 碰撞点
    double M_PI = 3.14159265358979323846;

    std::vector<igl::Hit> hits2;//碰撞检测
    std::vector<igl::Hit> hits3;
    int threebending = 0;//0顺时针极限1逆时针极限2正常折弯
    bool Is_Intersection[100] = { 0 };//避免碰撞部分----是否发生碰撞
    double Intersection_angle[100] = { 0.0 };//当前方向的limit
    double Intersection_angle2[100] = { 0.0 }; //相反方向的limit
    bool intersection_copy = 0;//检测碰撞
    bool angle_limit_change = 0;//
    bool threebending_init = 0;//初始化 只执行一次
    bool Intersection_al = 0;//其余部分会发生碰撞

    bool cal_angle = 0;//是否计算所需角度
    bool cal_angle_init = 0;//只进行一次
    Eigen::RowVector3d pos_copy;
    bool cal_angle_ok = 0;
    Eigen::RowVector3d b;//判断旋转角度是否完成

    bool is_rotating = false;//开始
    bool rotation_change = 0;//反方向执行一次
    Eigen::RowVector3d zeropoint;
    //std::vector<double> angles_close;

private:
    
};


class Collision_check {
public:
    Collision_check();
    
    bool checkCollision(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx);
    double Collision_check::angle_cal(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3);

    Eigen::RowVector3d wireStartPoint;
    const double M_PI = 3.14159265358979323846;

    int step = 0; //整体旋转步骤


    Eigen::RowVector3d center3;// 模型3中心
    Eigen::RowVector3d a;//判断旋转角度是否完成


    //int updownstep = 0;//上下步骤
    //int updownstepcount = 0;//上下多少
    //double z_offset = 0.5;//上下速度
    //double theta = 0.03;// 控制旋转角度和方向的变化
    //bool rotation_change = 0;//反方向执行一次
    //bool blenderdirchanged = 0;//是否未在范围内并且需要改变方向才能进入范围

    Eigen::MatrixXd V1, V2, V3;
    Eigen::MatrixXi F1, F2, F3;

    Eigen::RowVector3d min_corner;
    Eigen::RowVector3d max_corner;
    bool isIntersecting(const Eigen::Vector3d& point_start, const Eigen::Vector3d& point_end, const Eigen::RowVector3d& min_corner, const Eigen::RowVector3d& max_corner);

    std::vector<igl::Hit> hits;
    //Eigen::MatrixXi IF;
    //bool intersection = 0;//检测碰撞
    //Eigen::RowVector3d pos;//碰撞点
    //Eigen::RowVector3d pos2;//碰撞点之前位置
    //double intersection_1, intersection_2, intersection_3, intersection_4, intersection_5;  //ef 需要 double 碰撞点
    //std::vector<igl::Hit> hits2;//碰撞检测
    //std::vector<igl::Hit> hits3;
    //Eigen::RowVector3d V2point_center = (V2.row(112) + V2.row(205) + V2.row(173) + V2.row(148)) / 4;
    //

    //bool blenderback_place();//回撤完成位置
    //bool blenderback_place2();

    //double calculateAngle(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3);
    
    //bool blenderdirchange();//下一轮位置是否在完成位置范围内？大于小于需要改变方向否 结束后统一方向
    //Eigen::RowVector3d compute_rotation_point(const Eigen::RowVector3d& initial, const Eigen::RowVector3d& rotated, const Eigen::RowVector3d& initial2);//线旋转
    //void blenderback();//pin 下+回撤
    //void wireforward();//线走+pin 上
    //void wirebending();//弯曲
    //void drawpic();//绘制二维坐标轴
    ////bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);//控制输入
    //void drawAnimation();



    //
    //int anglesnum;
    //
    //bool blenderdirchanged = 0;//是否未在范围内并且需要改变方向才能进入范围

    //
    //
 
    //int already_length[100] = { 0 };
    


    
    //int threebending = 0;//0顺时针极限1逆时针极限2正常折弯
    //bool Is_Intersection[100] = { 0 };//避免碰撞部分----是否发生碰撞
    //double Intersection_angle[100] = { 0.0 };//当前方向的limit
    //double Intersection_angle2[100] = { 0.0 }; //相反方向的limit
    //bool intersection_copy = 0;//检测碰撞
    //bool angle_limit_change = 0;//
    //bool threebending_init = 0;//初始化 只执行一次
    //bool Intersection_al = 0;//其余部分会发生碰撞

    //bool cal_angle = 0;//是否计算所需角度
    //bool cal_angle_init = 0;//只进行一次
    //Eigen::RowVector3d pos_copy;
    //bool cal_angle_ok = 0;
    //Eigen::RowVector3d b;//判断旋转角度是否完成


    //bool rotation_change = 0;//反方向执行一次
    //Eigen::RowVector3d zeropoint;


private:

};


//bool checkCollision(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx);
std::pair<int, int> findLongestNonCollidingSegment(const std::vector<double>& lengths, const std::vector<double>& angles, int startIndex);
std::pair<int, int> findLongestNonCollidingSegment2(const std::vector<double>& lengths, const std::vector<double>& angles, int startIndex);

#endif





