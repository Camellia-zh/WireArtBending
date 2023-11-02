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
    bool blenderback_place();//�س����λ��
    bool blenderback_place2();
    bool blenderdirchange();//��һ��λ���Ƿ������λ�÷�Χ�ڣ�����С����Ҫ�ı䷽��� ������ͳһ����
    Eigen::RowVector3d compute_rotation_point(const Eigen::RowVector3d& initial, const Eigen::RowVector3d& rotated, const Eigen::RowVector3d& initial2);//����ת
    void blenderback();//pin ��+�س�
    void wireforward();//����+pin ��
    void wirebending();//����
    void drawpic();//���ƶ�ά������
    
    
    //bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);//��������
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
    Eigen::RowVector3d a;//�ж���ת�Ƕ��Ƿ����
    bool blenderdirchanged = 0;//�Ƿ�δ�ڷ�Χ�ڲ�����Ҫ�ı䷽����ܽ��뷶Χ
    int updownstep = 0;//���²���
    int updownstepcount = 0;//���¶���
    double z_offset = 0.5;//�����ٶ�
    double theta = 0.03;// ������ת�ǶȺͷ���ı仯
    Eigen::RowVector3d V2point_center = (V2.row(112) + V2.row(205) + V2.row(173) + V2.row(148)) / 4;
    Eigen::RowVector3d jichupoint;
    Eigen::RowVector3d center3;
    int already_length[100] = { 0 };
    std::vector<igl::Hit> hits;//�߶�
    Eigen::MatrixXi IF;
    bool intersection = 0;//�����ײ
    Eigen::RowVector3d pos;//��ײ��
    Eigen::RowVector3d pos2;//��ײ��֮ǰλ��
    double intersection_1, intersection_2, intersection_3, intersection_4, intersection_5;  //ef ��Ҫ double ��ײ��
    double M_PI = 3.14159265358979323846;

    std::vector<igl::Hit> hits2;//��ײ���
    std::vector<igl::Hit> hits3;
    int threebending = 0;//0˳ʱ�뼫��1��ʱ�뼫��2��������
    bool Is_Intersection[100] = { 0 };//������ײ����----�Ƿ�����ײ
    double Intersection_angle[100] = { 0.0 };//��ǰ�����limit
    double Intersection_angle2[100] = { 0.0 }; //�෴�����limit
    bool intersection_copy = 0;//�����ײ
    bool angle_limit_change = 0;//
    bool threebending_init = 0;//��ʼ�� ִֻ��һ��
    bool Intersection_al = 0;//���ಿ�ֻᷢ����ײ

    bool cal_angle = 0;//�Ƿ��������Ƕ�
    bool cal_angle_init = 0;//ֻ����һ��
    Eigen::RowVector3d pos_copy;
    bool cal_angle_ok = 0;
    Eigen::RowVector3d b;//�ж���ת�Ƕ��Ƿ����

    bool is_rotating = false;//��ʼ
    bool rotation_change = 0;//������ִ��һ��
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

    int step = 0; //������ת����


    Eigen::RowVector3d center3;// ģ��3����
    Eigen::RowVector3d a;//�ж���ת�Ƕ��Ƿ����


    //int updownstep = 0;//���²���
    //int updownstepcount = 0;//���¶���
    //double z_offset = 0.5;//�����ٶ�
    //double theta = 0.03;// ������ת�ǶȺͷ���ı仯
    //bool rotation_change = 0;//������ִ��һ��
    //bool blenderdirchanged = 0;//�Ƿ�δ�ڷ�Χ�ڲ�����Ҫ�ı䷽����ܽ��뷶Χ

    Eigen::MatrixXd V1, V2, V3;
    Eigen::MatrixXi F1, F2, F3;

    Eigen::RowVector3d min_corner;
    Eigen::RowVector3d max_corner;
    bool isIntersecting(const Eigen::Vector3d& point_start, const Eigen::Vector3d& point_end, const Eigen::RowVector3d& min_corner, const Eigen::RowVector3d& max_corner);

    std::vector<igl::Hit> hits;
    //Eigen::MatrixXi IF;
    //bool intersection = 0;//�����ײ
    //Eigen::RowVector3d pos;//��ײ��
    //Eigen::RowVector3d pos2;//��ײ��֮ǰλ��
    //double intersection_1, intersection_2, intersection_3, intersection_4, intersection_5;  //ef ��Ҫ double ��ײ��
    //std::vector<igl::Hit> hits2;//��ײ���
    //std::vector<igl::Hit> hits3;
    //Eigen::RowVector3d V2point_center = (V2.row(112) + V2.row(205) + V2.row(173) + V2.row(148)) / 4;
    //

    //bool blenderback_place();//�س����λ��
    //bool blenderback_place2();

    //double calculateAngle(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3);
    
    //bool blenderdirchange();//��һ��λ���Ƿ������λ�÷�Χ�ڣ�����С����Ҫ�ı䷽��� ������ͳһ����
    //Eigen::RowVector3d compute_rotation_point(const Eigen::RowVector3d& initial, const Eigen::RowVector3d& rotated, const Eigen::RowVector3d& initial2);//����ת
    //void blenderback();//pin ��+�س�
    //void wireforward();//����+pin ��
    //void wirebending();//����
    //void drawpic();//���ƶ�ά������
    ////bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier);//��������
    //void drawAnimation();



    //
    //int anglesnum;
    //
    //bool blenderdirchanged = 0;//�Ƿ�δ�ڷ�Χ�ڲ�����Ҫ�ı䷽����ܽ��뷶Χ

    //
    //
 
    //int already_length[100] = { 0 };
    


    
    //int threebending = 0;//0˳ʱ�뼫��1��ʱ�뼫��2��������
    //bool Is_Intersection[100] = { 0 };//������ײ����----�Ƿ�����ײ
    //double Intersection_angle[100] = { 0.0 };//��ǰ�����limit
    //double Intersection_angle2[100] = { 0.0 }; //�෴�����limit
    //bool intersection_copy = 0;//�����ײ
    //bool angle_limit_change = 0;//
    //bool threebending_init = 0;//��ʼ�� ִֻ��һ��
    //bool Intersection_al = 0;//���ಿ�ֻᷢ����ײ

    //bool cal_angle = 0;//�Ƿ��������Ƕ�
    //bool cal_angle_init = 0;//ֻ����һ��
    //Eigen::RowVector3d pos_copy;
    //bool cal_angle_ok = 0;
    //Eigen::RowVector3d b;//�ж���ת�Ƕ��Ƿ����


    //bool rotation_change = 0;//������ִ��һ��
    //Eigen::RowVector3d zeropoint;


private:

};


//bool checkCollision(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx);
std::pair<int, int> findLongestNonCollidingSegment(const std::vector<double>& lengths, const std::vector<double>& angles, int startIndex);
std::pair<int, int> findLongestNonCollidingSegment2(const std::vector<double>& lengths, const std::vector<double>& angles, int startIndex);

#endif





