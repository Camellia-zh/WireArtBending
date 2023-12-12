#include "WireBending.h"
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Intersections_3/Segment_3_Segment_3.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_3 Point_3;
typedef K::Segment_3 Segment_3;

int anglesNumber = 400;
extern std::chrono::duration<double> timeCC1;
extern std::chrono::duration<double> timeCC2;//cro碰撞检测时间

Wire::Wire() {
    // Initialize points matrix with 30 rows and 3 columns
    points.resize(anglesNumber, 3);
}

void Wire::initialize() {
    for (int i = 0; i < points.size(); i++) {
        points.row(i) << 498.2, -240.48, 185.0;
    }
}

void Wire::modifyPoint(int index, double x, double y, double z) {
    if (index >= 0 && index < points.rows()) {
        points.row(index) << x, y, z;
    }
    else {
        std::cout << "Invalid point index." << std::endl;
    }
}


Bending::Bending(Eigen::MatrixXd& V1_, Eigen::MatrixXd& V2_, Eigen::MatrixXd& V3_,
    Eigen::MatrixXi& F1_, Eigen::MatrixXi& F2_, Eigen::MatrixXi& F3_,
	Eigen::MatrixXd& wire_, std::vector<double>& distances_, std::vector<double>& angles_, int& anglesnum_, std::vector<double>& angles2_)
    : V1(V1_), V2(V2_), V3(V3_), F1(F1_), F2(F2_), F3(F3_),
    w(wire_), distances(distances_), angles(angles_), anglesnum(anglesnum_), angles2(angles2_),
	a(30, 0, 0), jichupoint(498.2, -240.48, 185.0), center3(492.2, -240.46, 168.8),
	zeropoint(600, -200, 185.0) {}


double Bending::calculateAngle(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3) {
	// Calculate vectors between points
	Eigen::Vector3d vector1 = p1 - p2;
	Eigen::Vector3d vector2 = p3 - p2;
	// Calculate dot product of the two vectors
	double dotProduct = vector1.dot(vector2);
	// Calculate magnitudes of the vectors
	double magnitude1 = vector1.norm();
	double magnitude2 = vector2.norm();
	// Calculate the cosine of the angle between the vectors
	double cosine = dotProduct / (magnitude1 * magnitude2);
	// Calculate the angle in radians
	double angleInRadians = std::acos(cosine);
	// Convert angle from radians to degrees
	double angleInDegrees = angleInRadians * 180.0 / M_PI;

	return angleInDegrees;
}

bool Bending::blenderback_place() {
	double v1_x = V2point_center.x() - 483.2;
	double v1_y = V2point_center.y() - (-240.48);
	double v2_x = a.x();
	double v2_y = a.y();
	double cross = v1_x * v2_y - v1_y * v2_x;
	double dot = v1_x * v2_x + v1_y * v2_y;
	double norm_v1 = std::sqrt(v1_x * v1_x + v1_y * v1_y);
	double norm_v2 = std::sqrt(v2_x * v2_x + v2_y * v2_y);
	double angle = round(std::acos(dot / (norm_v1 * norm_v2)) * 180 / 3.1415926535);
	if (angles[step / 3] > 0)
	{	
		if (cross < 0 && ((angle > 60) && (angle < 100)))
		{
			return true;
		}
	}
	else if (angles[step / 3] < 0) {
		if (cross > 0 && ((angle > 80) && (angle < 100)))
		{
			return true;
		}
	}
	else
	{
		viewer.core().is_animating = !viewer.core().is_animating;
		return false;
	}
	return false;
}

bool Bending::blenderback_place2() {
	double v1_x = V2point_center.x() - 483.2;
	double v1_y = V2point_center.y() - (-240.48);
	double v2_x = a.x();
	double v2_y = a.y();
	double cross = v1_x * v2_y - v1_y * v2_x;
	double dot = v1_x * v2_x + v1_y * v2_y;
	double norm_v1 = std::sqrt(v1_x * v1_x + v1_y * v1_y);
	double norm_v2 = std::sqrt(v2_x * v2_x + v2_y * v2_y);
	double angle = round(std::acos(dot / (norm_v1 * norm_v2)) * 180 / 3.1415926535);
	if (angles[step / 3] < 0)
	{
		if (cross < 0 && ((angle > 60) && (angle < 100)))
		{
			return true;
		}
	}
	else if (angles[step / 3] > 0) {
		if (cross > 0 && ((angle > 80) && (angle < 100)))
		{
			return true;
		}
	}
	else
	{
		viewer.core().is_animating = !viewer.core().is_animating;
		return false;
	}
	return false;
}

bool Bending::blenderdirchange() {
	double v1_x = V2point_center.x() - 483.2;
	double v1_y = V2point_center.y() - (-240.48);
	double v2_x = a.x();
	double v2_y = a.y();
	double cross = v1_x * v2_y - v1_y * v2_x;
	double dot = v1_x * v2_x + v1_y * v2_y;
	double norm_v1 = std::sqrt(v1_x * v1_x + v1_y * v1_y);
	double norm_v2 = std::sqrt(v2_x * v2_x + v2_y * v2_y);
	double angle = round(std::acos(dot / (norm_v1 * norm_v2)) * 180 / 3.1415926535);
	if (angle > 100)
	{
		blenderdirchanged = 1;
		return true;
	}
	return false;
}

Eigen::RowVector3d Bending::compute_rotation_point(const Eigen::RowVector3d& initial, const Eigen::RowVector3d& rotated, const Eigen::RowVector3d& initial2) {

	//Eigen::RowVector3d initial_rel = initial - jichupoint;  // 计算initial相对于jichupoint的位置
	//Eigen::RowVector3d rotated_rel = rotated - jichupoint;  // 计算rotated相对于jichupoint的位置
	//Eigen::RowVector3d axis = initial_rel.cross(rotated_rel);  // 计算旋转轴
	//axis.normalize();  // 单位化旋转轴
	//double angle = acos(initial_rel.dot(rotated_rel) / (initial_rel.norm() * rotated_rel.norm()));  // 计算旋转角度
	//// 用Rodrigues' rotation formula计算旋转矩阵
	//Eigen::Matrix3d K;
	//K << 0, -axis(2), axis(1),
	//	axis(2), 0, -axis(0),
	//	-axis(1), axis(0), 0;
	//Eigen::Matrix3d rotationMatrix = Eigen::Matrix3d::Identity() + sin(angle) * K + (1 - cos(angle)) * K * K;
	//Eigen::RowVector3d point_rel = initial2 - jichupoint;  // 计算point相对于jichupoint的位置
	//Eigen::RowVector3d point_rotated_rel = rotationMatrix * point_rel.transpose();  // 旋转point
	//Eigen::RowVector3d point_rotated = point_rotated_rel + jichupoint;  // 从相对位置回到绝对位置
	//return point_rotated;

	Eigen::Vector3d p1_initial_3d = (initial - jichupoint).transpose();
	Eigen::Vector3d p1_rotated_3d = (rotated - jichupoint).transpose();
	Eigen::Vector3d p2_initial_3d = (initial2 - jichupoint).transpose();
	Eigen::Vector2d p1_initial(p1_initial_3d.x(), p1_initial_3d.y());
	Eigen::Vector2d p1_rotated(p1_rotated_3d.x(), p1_rotated_3d.y());
	Eigen::Vector2d p2_initial(p2_initial_3d.x(), p2_initial_3d.y());
	double dot = p1_initial.dot(p1_rotated);
	double det = p1_initial.x() * p1_rotated.y() - p1_initial.y() * p1_rotated.x();
	double angle = atan2(det, dot);
	Eigen::Matrix2d rotationMatrix;
	rotationMatrix << cos(angle), -sin(angle),
		sin(angle), cos(angle);
	p2_initial = rotationMatrix * p2_initial;
	Eigen::Vector3d p2_rotated_3d(p2_initial.x(), p2_initial.y(), p2_initial_3d.z());
	return (p2_rotated_3d).transpose() + jichupoint;

	//// 计算旋转轴
	//Eigen::Vector3d rotation_axis = p1_initial.cross(p1_rotated).normalized();
	//// 计算旋转角度
	//double cos_theta = p1_initial.normalized().dot(p1_rotated.normalized());
	//double theta = acos(cos_theta);
	//// 计算旋转矩阵
	//Eigen::Matrix3d K = Eigen::Matrix3d::Zero();
	//K(0, 1) = -rotation_axis(2);
	//K(0, 2) = rotation_axis(1);
	//K(1, 0) = rotation_axis(2);
	//K(1, 2) = -rotation_axis(0);
	//K(2, 0) = -rotation_axis(1);
	//K(2, 1) = rotation_axis(0);
	//Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + sin(theta) * K + (1 - cos_theta) * K * K;
	//return (R * p2_initial).transpose() + jichupoint;
}

void Bending::blenderback() {
	if (updownstep == 0)
	{
		z_offset = -0.5;
		updownstepcount++;
		V2 = V2.rowwise() + Eigen::Vector3d(0.0, 0.0, z_offset).transpose();
		viewer.data_list[1].set_vertices(V2);
		if (updownstepcount == 50)
		{
			if (step == 0)
			{
				updownstep += 1;
			}
			else {
				if (step / 3 < anglesnum)//没到最后一个弯
				{
					updownstep += 1;
					if ((angles[step / 3] * angles[(step / 3) - 1]) > 0)
					{
						theta = -theta;
					}
					else if ((angles[step / 3] * angles[(step / 3) - 1]) < 0)
					{
						V2point_center = (V2.row(112) + V2.row(205) + V2.row(173) + V2.row(148)) / 4;
						if (blenderdirchange())
						{
							theta = -theta;
						}
					}
					else
					{
						viewer.core().is_animating = !viewer.core().is_animating;
					}
				}
				else {
					step++;
				}
			}
		}
	}
	if (updownstep == 1)
	{
		// Update the rotation angle
		Eigen::Matrix3d rotation;
		rotation = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
		// Apply the rotation to the rotating model
		V2 = (V2.rowwise() - center3) * rotation;
		V2.rowwise() += center3;
		V3 = (V3.rowwise() - center3) * rotation;
		V3.rowwise() += center3;
		// Update the viewer with the rotated model
		viewer.data_list[1].set_vertices(V2);
		viewer.data_list[2].set_vertices(V3);
		V2point_center = (V2.row(112) + V2.row(205) + V2.row(173) + V2.row(148)) / 4;
		if (blenderback_place())
		{
			updownstep = 0;
			updownstepcount = 0;//切换下一步时，updownstep updownstepcount为0
			step++;
			if (blenderdirchanged)
			{
				theta = -theta;
				blenderdirchanged = 0;
			}
		}
	}
}

void Bending::wireforward() {
	if (updownstep == 0)
	{
		for (int i = 0; i <= step / 3; i++)
		{
			w(i, 0) += 1;
		}
		already_length[step / 3]++;
		viewer.data().clear_edges();
		for (int j = 0; j <= step / 3; j++)
		{
			viewer.data().add_edges(w.row(j), w.row(j + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
		}
		if (already_length[step / 3] >= distances[step / 3])
		{
			updownstep++;
		}
	}
	if (updownstep == 1)
	{
		if (step / 3 >= anglesnum)
		{
			viewer.core().is_animating = !viewer.core().is_animating;
		}
		z_offset = +0.5;
		updownstepcount++;
		V2 = V2.rowwise() + Eigen::Vector3d(0.0, 0.0, z_offset).transpose();
		viewer.data_list[1].set_vertices(V2);
		if (updownstepcount == 50)
		{
			step++;
			updownstep = 0;
			updownstepcount = 0;
			theta = -theta;
			if (angles[step / 3] > 0 && step == 2)//step 已经加加了
			{
				theta = -theta;
			}
		}
	}
}

void Bending::wirebending() {
	// Update the rotation angle
	Eigen::Matrix3d rotation;
	rotation = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
	// Apply the rotation to the rotating model
	V2 = (V2.rowwise() - center3) * rotation;
	V2.rowwise() += center3;
	V3 = (V3.rowwise() - center3) * rotation;
	V3.rowwise() += center3;
	// Update the viewer with the rotated model
	viewer.data_list[1].set_vertices(V2);
	viewer.data_list[2].set_vertices(V3);
	// Compute the intersection between the ray and the model
	if (igl::ray_mesh_intersect(w.row((step / 3) + 1), w.row(step / 3) - w.row((step / 3) + 1), V2, F2, hits))
	{
		if (!intersection)
		{
			intersection = 1;
			w.row(step / 3) = V2.row(F2(hits[0].id, 0)) * (1 - hits[0].u - hits[0].v) + V2.row(F2(hits[0].id, 1)) * hits[0].u + V2.row(F2(hits[0].id, 2)) * hits[0].v;
			intersection_1 = F2(hits[0].id, 0);
			intersection_2 = F2(hits[0].id, 1);
			intersection_3 = F2(hits[0].id, 2);
			intersection_4 = hits[0].u;
			intersection_5 = hits[0].v;
			pos2 = w.row(step / 3);
			//cout << "intersection_1: " << intersection_1 << endl;
		}
		pos = V2.row(intersection_1) * (1 - intersection_4 - intersection_5) + V2.row(intersection_2) * intersection_4 + V2.row(intersection_3) * intersection_5;
		w.row(step / 3) = pos;
		w(step / 3, 0) = w((step / 3) + 1, 0) + (distances[step / 3] * (pos.x() - w((step / 3) + 1, 0)) / sqrt(pow(pos.x() - w((step / 3) + 1, 0), 2) + pow(pos.y() - w((step / 3) + 1, 1), 2)));
		w(step / 3, 1) = w((step / 3) + 1, 1) + (distances[step / 3] * (pos.y() - w((step / 3) + 1, 1)) / sqrt(pow(pos.x() - w((step / 3) + 1, 0), 2) + pow(pos.y() - w((step / 3) + 1, 1), 2)));
		//cout << w.row(step / 3) << endl;
		if (intersection)
		{
			for (int j = 0; j < step / 3; j++)
			{
				Eigen::RowVector3d p2_rotated = compute_rotation_point(pos2, pos, w.row(j));
				w.row(j) = p2_rotated;
			}
			pos2 = pos;
		}
	}
	viewer.data().clear_edges();
	for (int i = 0; i < step / 3; i++)
	{
		bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits2) && ((hits2[0].t < 1) && (hits2[0].t > 0));
		bool wireinter2 = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V2, F2, hits3) && ((hits3[0].t < 1) && (hits3[0].t > 0));
		if (wireinter || wireinter2)
		{
			viewer.data().add_edges(w.row(i), w.row(i + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
		}
		else {
			viewer.data().add_edges(w.row(i), w.row(i + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
		}
		if (wireinter && (!Intersection_al))
		{
			Is_Intersection[step / 3] = 1;
			Eigen::RowVector3d b2 = w.row(step / 3) - w.row((step / 3) + 1);
			double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
			double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
			double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
			double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
			double angle2 = acos(cos_angle2); // 夹角的弧度值
			double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
			if (angles[step / 3] > 0)
			{
				Intersection_angle[step / 3] = angle_degrees2;
			}
			else {
				Intersection_angle[step / 3] = -angle_degrees2;
			}
			Intersection_al = 1;
		}
	}
	//viewer.data().add_edges(w.row(step / 3), w.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
	// Check for collisions
	bool collision = igl::copyleft::cgal::intersect_other(V1, F1, V2, F2, false, IF);
	// If a collision is detected, reverse the rotation direction
	if (collision) {
		theta = -theta;
	}
	//
	b = w.row(step / 3) - w.row((step / 3) + 1);
	double dot_product = a.x() * b.x() + a.y() * b.y(); // 向量点乘
	double a_length = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
	double b_length = sqrt(b.x() * b.x() + b.y() * b.y()); // 向量 cd 的长度
	double cos_angle = dot_product / (a_length * b_length); // 夹角的余弦值
	double angle = acos(cos_angle); // 夹角的弧度值
	double angle_degrees = angle * 180 / M_PI; // 将弧度值转换为角度值
	if (angle_degrees >= 130)
	{
		Is_Intersection[step / 3] = 1;
		viewer.data().add_edges(w.row(step / 3), w.row((step / 3) + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
	}
	else
	{
		viewer.data().add_edges(w.row(step / 3), w.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
	}
	if (angle_degrees >= abs(angles[step / 3]))
	{
		step++;
		intersection = 0;
		Intersection_al = 0;
		threebending = 0;
	}
}

void Bending::drawpic() {
	Eigen::RowVector3d yuandian = zeropoint;
	//viewer.data().add_points(yuandian, Eigen::RowVector3d(1, 1, 1));
	viewer.data().add_edges(yuandian, yuandian + Eigen::RowVector3d(480, 0, 0), Eigen::RowVector3d(1, 1, 1));
	viewer.data().add_edges(yuandian + Eigen::RowVector3d(480, 0, 0), yuandian + Eigen::RowVector3d(480, 0, 0) + Eigen::RowVector3d(-3, 0, 3), Eigen::RowVector3d(1, 1, 1));
	viewer.data().add_edges(yuandian + Eigen::RowVector3d(480, 0, 0), yuandian + Eigen::RowVector3d(480, 0, 0) + Eigen::RowVector3d(-3, 0, -3), Eigen::RowVector3d(1, 1, 1));
	viewer.data().add_edges(yuandian - Eigen::RowVector3d(0, 0, 180), yuandian + Eigen::RowVector3d(0, 0, 180), Eigen::RowVector3d(1, 1, 1));
	viewer.data().add_edges(yuandian + Eigen::RowVector3d(0, 0, 180), yuandian + Eigen::RowVector3d(0, 0, 180) + Eigen::RowVector3d(-3, 0, -3), Eigen::RowVector3d(1, 1, 1));
	viewer.data().add_edges(yuandian + Eigen::RowVector3d(0, 0, 180), yuandian + Eigen::RowVector3d(0, 0, 180) + Eigen::RowVector3d(3, 0, -3), Eigen::RowVector3d(1, 1, 1));
	int al_length_head = 0;
	int al_length_rear = 0;
	for (int i = 0; i < step / 3; i++)
	{
		al_length_rear += distances[i];
		Eigen::RowVector3d shiji = yuandian + Eigen::RowVector3d(al_length_rear, 0, angles[i]);
		if (Is_Intersection[i])
		{
			viewer.data().add_points(shiji, Eigen::RowVector3d(1, 0, 0));//红
			angles2[i] = Intersection_angle[i]-2;
			viewer.data().add_edges(yuandian + Eigen::RowVector3d(al_length_head, 0, (Intersection_angle[i])), yuandian + Eigen::RowVector3d(al_length_rear, 0, (Intersection_angle[i])), Eigen::RowVector3d(0, 1, 0));
			viewer.data().add_edges(yuandian + Eigen::RowVector3d(al_length_head, 0, (Intersection_angle2[i])), yuandian + Eigen::RowVector3d(al_length_rear, 0, (Intersection_angle2[i])), Eigen::RowVector3d(0, 1, 0));
		}
		else {
			viewer.data().add_points(shiji, Eigen::RowVector3d(0, 1, 0));//绿
			angles2[i] = angles[i];
			viewer.data().add_edges(yuandian + Eigen::RowVector3d(al_length_head, 0, (Intersection_angle[i])), yuandian + Eigen::RowVector3d(al_length_rear, 0, (Intersection_angle[i])), Eigen::RowVector3d(0, 1, 0));
			viewer.data().add_edges(yuandian + Eigen::RowVector3d(al_length_head, 0, (Intersection_angle2[i])), yuandian + Eigen::RowVector3d(al_length_rear, 0, (Intersection_angle2[i])), Eigen::RowVector3d(0, 1, 0));
		}
		al_length_head += distances[i];
		/*cout << "Intersection_angle[i]:" << Intersection_angle[i] << " angles[i]:" << angles[i] << endl;
		viewer.core().is_animating = !viewer.core().is_animating;*/

		/*std::stringstream l1;
		l1 << zeropoint(0) << ", " << zeropoint(1) << ", " << zeropoint(2);
		viewer.data().add_label(zeropoint + Eigen::RowVector3d(10, 10, 0), l1.str());
		viewer.data().show_custom_labels = true;*/
	}
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
	switch (key)
	{
	case ' ':
		viewer.core().is_animating = !viewer.core().is_animating;
		return true;
	}
	return false;
}

void Bending::drawAnimation() {
	// 模型1静态
	viewer.data().set_mesh(V1, F1);
	viewer.data().set_face_based(true);  //设置网格的着色模式 基于面的着色模式
	viewer.data().show_lines = false;
	// 模型23旋转
	viewer.append_mesh();
	viewer.data().set_mesh(V2, F2);
	viewer.data().set_face_based(true);
	viewer.data().show_lines = false;
	viewer.append_mesh();
	viewer.data().set_mesh(V3, F3);
	viewer.data().set_face_based(true);
	viewer.data().show_lines = false;
	// 使用add_edges函数添加线段，其中第一个参数表示起点坐标，第二个参数表示终点坐标，第三个参数表示线段颜色
	viewer.data().add_edges(w.row(0), w.row(1), Eigen::RowVector3d(0.0, 1.0, 0.0));

	Eigen::MatrixXd V2_copy = V2;
	Eigen::MatrixXi F2_copy = F2;
	viewer.append_mesh();
	viewer.data().set_mesh(V2_copy, F2_copy);
	viewer.data().show_lines = false;
	Eigen::MatrixXd V3_copy = V3;
	Eigen::MatrixXi F3_copy = F3;
	viewer.append_mesh();
	viewer.data().set_mesh(V3_copy, F3_copy);
	viewer.data().show_lines = false;
	viewer.data_list[3].show_faces = false;
	viewer.data_list[4].show_faces = false;
	Eigen::MatrixXd w_copy = w;
	double theta_copy = theta;

	double point_size = 5.0;
	viewer.data().point_size = point_size;
	double line_width = 3.0;
	viewer.data().line_width = line_width;

	// 预绘制回调以更新旋转模型并检查碰撞
	viewer.callback_pre_draw = [&](igl::opengl::glfw::Viewer& viewer)-> bool
	{
		// Check if the rotation is enabled
		if (viewer.core().is_animating)
		{
			if (step % 3 == 0)
			{
				blenderback();
			}
			if (step % 3 == 1)
			{
				wireforward();
			}
			if (step % 3 == 2)
			{
				if (threebending == 0)
				{
					if (!threebending_init)//只进行一次
					{
						viewer.data_list[1].show_faces = false;
						viewer.data_list[2].show_faces = false;
						threebending_init = true;
						V2_copy = V2; F2_copy = F2; V3_copy = V3; F3_copy = F3; w_copy = w; theta_copy = theta;
						viewer.data_list[3].show_faces = true;
						viewer.data_list[4].show_faces = true;
					}
					// Update the rotation angle
					Eigen::Matrix3d rotation;
					rotation = Eigen::AngleAxisd(theta_copy, Eigen::Vector3d::UnitZ());
					// Apply the rotation to the rotating model
					V2_copy = (V2_copy.rowwise() - center3) * rotation;
					V2_copy.rowwise() += center3;
					V3_copy = (V3_copy.rowwise() - center3) * rotation;
					V3_copy.rowwise() += center3;
					// Update the viewer with the rotated model
					viewer.data_list[3].set_vertices(V2_copy);
					viewer.data_list[4].set_vertices(V3_copy);
					// Compute the intersection between the ray and the model
					if (igl::ray_mesh_intersect(w_copy.row((step / 3) + 1), w_copy.row(step / 3) - w_copy.row((step / 3) + 1), V2_copy, F2_copy, hits))
					{
						if (!intersection_copy)
						{
							intersection_copy = 1;
							w_copy.row(step / 3) = V2_copy.row(F2_copy(hits[0].id, 0)) * (1 - hits[0].u - hits[0].v) + V2_copy.row(F2_copy(hits[0].id, 1)) * hits[0].u + V2_copy.row(F2_copy(hits[0].id, 2)) * hits[0].v;
							intersection_1 = F2_copy(hits[0].id, 0);
							intersection_2 = F2_copy(hits[0].id, 1);
							intersection_3 = F2_copy(hits[0].id, 2);
							intersection_4 = hits[0].u;
							intersection_5 = hits[0].v;
							pos2 = w_copy.row(step / 3);
						}
						pos = V2_copy.row(intersection_1) * (1 - intersection_4 - intersection_5) + V2_copy.row(intersection_2) * intersection_4 + V2_copy.row(intersection_3) * intersection_5;
						w_copy.row(step / 3) = pos;
						w_copy(step / 3, 0) = w_copy((step / 3) + 1, 0) + (distances[step / 3] * (pos.x() - w_copy((step / 3) + 1, 0)) / sqrt(pow(pos.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos.y() - w_copy((step / 3) + 1, 1), 2)));
						w_copy(step / 3, 1) = w_copy((step / 3) + 1, 1) + (distances[step / 3] * (pos.y() - w_copy((step / 3) + 1, 1)) / sqrt(pow(pos.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos.y() - w_copy((step / 3) + 1, 1), 2)));
						if (intersection_copy)
						{
							for (int j = 0; j < step / 3; j++)
							{
								Eigen::RowVector3d p2_rotated = compute_rotation_point(pos2, pos, w_copy.row(j));
								w_copy.row(j) = p2_rotated;
							}
							pos2 = pos;
						}
					}
					viewer.data().clear_edges();
					for (int i = 0; i < step / 3; i++)
					{
						bool wireinter = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V1, F1, hits2) && ((hits2[0].t < 1) && (hits2[0].t > 0));
						bool wireinter2 = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V2_copy, F2_copy, hits3) && ((hits3[0].t < 1) && (hits3[0].t > 0));
						if (wireinter || wireinter2)
						{
							//Is_Intersection[step / 3] = 1;
							Eigen::RowVector3d b2 = w_copy.row(step / 3) - w_copy.row((step / 3) + 1);
							double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
							double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
							double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
							double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
							double angle2 = acos(cos_angle2); // 夹角的弧度值
							double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
							if (angles[step / 3] > 0)
							{
								Intersection_angle[step / 3] = angle_degrees2;
								angle_limit_change = 1;
							}
							else {
								Intersection_angle[step / 3] = -angle_degrees2;
								angle_limit_change = 1;
							}
							viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
						}
						else {
							viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
						}
					}
					Eigen::RowVector3d b2 = w_copy.row(step / 3) - w_copy.row((step / 3) + 1);
					double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
					double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
					double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
					double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
					double angle2 = acos(cos_angle2); // 夹角的弧度值
					double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
					if (angle_degrees2 >= 130)
					{
						viewer.data().add_edges(w_copy.row(step / 3), w_copy.row((step / 3) + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
						if (angles[step / 3] > 0)
						{
							//Is_Intersection[step / 3] = 1;
							Intersection_angle[step / 3] = angle_degrees2;
							angle_limit_change = 1;
						}
						else {
							//Is_Intersection[step / 3] = 1;
							Intersection_angle[step / 3] = -angle_degrees2;
							angle_limit_change = 1;
						}
					}
					else
					{
						viewer.data().add_edges(w_copy.row(step / 3), w_copy.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
					}
					if (angle_limit_change)
					{
						threebending += 1;
						threebending_init = 0;
						intersection_copy = 0;
						angle_limit_change = 0;
					}
				}
				else if (threebending == 1)
				{
					if (!threebending_init)//只进行一次
					{
						threebending_init = true;
						V2_copy = V2; F2_copy = F2; V3_copy = V3; F3_copy = F3; w_copy = w; theta_copy = theta;
						viewer.data_list[3].show_faces = true;
						viewer.data_list[4].show_faces = true;
						viewer.data_list[3].set_vertices(V2_copy);
						viewer.data_list[4].set_vertices(V3_copy);
						viewer.data().clear_edges();
					}
					for (int j = 0; j <= step / 3; j++)
					{
						viewer.data().add_edges(w.row(j), w.row(j + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
					}
					if (updownstep == 0)
					{
						z_offset = -0.5;
						updownstepcount++;
						V2_copy = V2_copy.rowwise() + Eigen::Vector3d(0.0, 0.0, z_offset).transpose();
						viewer.data_list[3].set_vertices(V2_copy);
						if (updownstepcount == 50)
						{
							updownstep += 1;
						}
					}
					if (updownstep == 1)
					{
						// Update the rotation angle
						Eigen::Matrix3d rotation;
						rotation = Eigen::AngleAxisd(theta_copy, Eigen::Vector3d::UnitZ());
						// Apply the rotation to the rotating model
						V2_copy = (V2_copy.rowwise() - center3) * rotation;
						V2_copy.rowwise() += center3;
						V3_copy = (V3_copy.rowwise() - center3) * rotation;
						V3_copy.rowwise() += center3;
						// Update the viewer with the rotated model
						viewer.data_list[3].set_vertices(V2_copy);
						viewer.data_list[4].set_vertices(V3_copy);
						V2point_center = (V2_copy.row(112) + V2_copy.row(205) + V2_copy.row(173) + V2_copy.row(148)) / 4;
						if (blenderback_place2())
						{
							updownstep++;
							updownstepcount = 0;//切换下一步时，updownstep updownstepcount为0
							theta_copy = -theta_copy;
						}
					}
					if (updownstep == 2)
					{
						if (step / 3 >= anglesnum)
						{
							viewer.core().is_animating = !viewer.core().is_animating;
						}
						z_offset = +0.5;
						updownstepcount++;
						V2_copy = V2_copy.rowwise() + Eigen::Vector3d(0.0, 0.0, z_offset).transpose();
						viewer.data_list[3].set_vertices(V2_copy);
						if (updownstepcount == 50)
						{
							updownstep++;
							updownstepcount = 0;
						}
					}
					if (updownstep == 3)
					{
						// Update the rotation angle
						Eigen::Matrix3d rotation;
						rotation = Eigen::AngleAxisd(theta_copy, Eigen::Vector3d::UnitZ());
						// Apply the rotation to the rotating model
						V2_copy = (V2_copy.rowwise() - center3) * rotation;
						V2_copy.rowwise() += center3;
						V3_copy = (V3_copy.rowwise() - center3) * rotation;
						V3_copy.rowwise() += center3;
						// Update the viewer with the rotated model
						viewer.data_list[3].set_vertices(V2_copy);
						viewer.data_list[4].set_vertices(V3_copy);
						// Compute the intersection between the ray and the model
						if (igl::ray_mesh_intersect(w_copy.row((step / 3) + 1), w_copy.row(step / 3) - w_copy.row((step / 3) + 1), V2_copy, F2_copy, hits))
						{
							if (!intersection_copy)
							{
								intersection_copy = 1;
								w_copy.row(step / 3) = V2_copy.row(F2_copy(hits[0].id, 0)) * (1 - hits[0].u - hits[0].v) + V2_copy.row(F2_copy(hits[0].id, 1)) * hits[0].u + V2_copy.row(F2_copy(hits[0].id, 2)) * hits[0].v;
								intersection_1 = F2_copy(hits[0].id, 0);
								intersection_2 = F2_copy(hits[0].id, 1);
								intersection_3 = F2_copy(hits[0].id, 2);
								intersection_4 = hits[0].u;
								intersection_5 = hits[0].v;
								pos2 = w_copy.row(step / 3);
							}
							pos = V2_copy.row(intersection_1) * (1 - intersection_4 - intersection_5) + V2_copy.row(intersection_2) * intersection_4 + V2_copy.row(intersection_3) * intersection_5;
							w_copy.row(step / 3) = pos;
							w_copy(step / 3, 0) = w_copy((step / 3) + 1, 0) + (distances[step / 3] * (pos.x() - w_copy((step / 3) + 1, 0)) / sqrt(pow(pos.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos.y() - w_copy((step / 3) + 1, 1), 2)));
							w_copy(step / 3, 1) = w_copy((step / 3) + 1, 1) + (distances[step / 3] * (pos.y() - w_copy((step / 3) + 1, 1)) / sqrt(pow(pos.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos.y() - w_copy((step / 3) + 1, 1), 2)));
							if (intersection_copy)
							{
								for (int j = 0; j < step / 3; j++)
								{
									Eigen::RowVector3d p2_rotated = compute_rotation_point(pos2, pos, w_copy.row(j));
									w_copy.row(j) = p2_rotated;
								}
								pos2 = pos;
							}
						}
						viewer.data().clear_edges();
						for (int i = 0; i < step / 3; i++)
						{
							bool wireinter = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V1, F1, hits2) && ((hits2[0].t < 1) && (hits2[0].t > 0));
							bool wireinter2 = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V2_copy, F2_copy, hits3) && ((hits3[0].t < 1) && (hits3[0].t > 0));
							if (wireinter || wireinter2)
							{
								Eigen::RowVector3d b2 = w_copy.row(step / 3) - w_copy.row((step / 3) + 1);
								double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
								double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
								double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
								double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
								double angle2 = acos(cos_angle2); // 夹角的弧度值
								double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
								if (angles[step / 3] > 0)
								{
									Intersection_angle[step / 3] = -angle_degrees2;
									angle_limit_change = 1;
								}
								else {
									Intersection_angle2[step / 3] = angle_degrees2;
									angle_limit_change = 1;
								}
								viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
							}
							else {
								viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
							}
						}
						Eigen::RowVector3d b2 = w_copy.row(step / 3) - w_copy.row((step / 3) + 1);
						double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
						double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
						double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
						double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
						double angle2 = acos(cos_angle2); // 夹角的弧度值
						double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
						if (angle_degrees2 >= 130)
						{
							viewer.data().add_edges(w_copy.row(step / 3), w_copy.row((step / 3) + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
							if (angles[step / 3] > 0)
							{
								Intersection_angle2[step / 3] = -angle_degrees2;
								angle_limit_change = 1;
							}
							else {
								Intersection_angle2[step / 3] = angle_degrees2;
								angle_limit_change = 1;

							}
						}
						else
						{
							viewer.data().add_edges(w_copy.row(step / 3), w_copy.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
						}
						if (angle_limit_change)
						{
							threebending++;
							threebending_init = 0;
							intersection_copy = 0;
							angle_limit_change = 0;
							updownstep = 0;
							viewer.data_list[1].show_faces = true;
							viewer.data_list[2].show_faces = true;
							viewer.data_list[3].show_faces = false;
							viewer.data_list[4].show_faces = false;
						}
					}
				}
				else if (threebending == 2)
				{
					// Update the rotation angle
					Eigen::Matrix3d rotation;
					rotation = Eigen::AngleAxisd(theta, Eigen::Vector3d::UnitZ());
					// Apply the rotation to the rotating model
					V2 = (V2.rowwise() - center3) * rotation;
					V2.rowwise() += center3;
					V3 = (V3.rowwise() - center3) * rotation;
					V3.rowwise() += center3;
					// Update the viewer with the rotated model
					viewer.data_list[1].set_vertices(V2);
					viewer.data_list[2].set_vertices(V3);
					// Compute the intersection between the ray and the model
					if (igl::ray_mesh_intersect(w.row((step / 3) + 1), w.row(step / 3) - w.row((step / 3) + 1), V2, F2, hits))
					{
						if (!intersection)
						{
							intersection = 1;
							w.row(step / 3) = V2.row(F2(hits[0].id, 0)) * (1 - hits[0].u - hits[0].v) + V2.row(F2(hits[0].id, 1)) * hits[0].u + V2.row(F2(hits[0].id, 2)) * hits[0].v;
							intersection_1 = F2(hits[0].id, 0);
							intersection_2 = F2(hits[0].id, 1);
							intersection_3 = F2(hits[0].id, 2);
							intersection_4 = hits[0].u;
							intersection_5 = hits[0].v;
							pos2 = w.row(step / 3);
							//cout << "intersection_1: " << intersection_1 << endl;
							//if (cal_angle)
							//{
							//	V2_copy = V2; F2_copy = F2; V3_copy = V3; F3_copy = F3; w_copy = w; theta_copy = theta;
							//	//cout << w_copy(step / 3, 0) <<"  "<< w_copy(step / 3, 1) << endl;
							//	w_copy(step / 3, 0) = w_copy((step / 3) + 1, 0) + (distances[step / 3] * (pos2.x() - w_copy((step / 3) + 1, 0)) / sqrt(pow(pos2.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos2.y() - w_copy((step / 3) + 1, 1), 2)));
							//	w_copy(step / 3, 1) = w_copy((step / 3) + 1, 1) + (distances[step / 3] * (pos2.y() - w_copy((step / 3) + 1, 1)) / sqrt(pow(pos2.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos2.y() - w_copy((step / 3) + 1, 1), 2)));
							//}
						}
						pos = V2.row(intersection_1) * (1 - intersection_4 - intersection_5) + V2.row(intersection_2) * intersection_4 + V2.row(intersection_3) * intersection_5;
						w.row(step / 3) = pos;
						w(step / 3, 0) = w((step / 3) + 1, 0) + (distances[step / 3] * (pos.x() - w((step / 3) + 1, 0)) / sqrt(pow(pos.x() - w((step / 3) + 1, 0), 2) + pow(pos.y() - w((step / 3) + 1, 1), 2)));
						w(step / 3, 1) = w((step / 3) + 1, 1) + (distances[step / 3] * (pos.y() - w((step / 3) + 1, 1)) / sqrt(pow(pos.x() - w((step / 3) + 1, 0), 2) + pow(pos.y() - w((step / 3) + 1, 1), 2)));
						//cout << w.row(step / 3) << endl;
						if (intersection)
						{
							for (int j = 0; j < step / 3; j++)
							{
								Eigen::RowVector3d p2_rotated = compute_rotation_point(pos2, pos, w.row(j));
								w.row(j) = p2_rotated;
							}
							pos2 = pos;
						}
					}
					viewer.data().clear_edges();
					for (int i = 0; i < step / 3; i++)
					{
						bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits2) && ((hits2[0].t < 1) && (hits2[0].t > 0));
						bool wireinter2 = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V2, F2, hits3) && ((hits3[0].t < 1) && (hits3[0].t > 0));
						if (wireinter || wireinter2)
						{
							viewer.data().add_edges(w.row(i), w.row(i + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
						}
						else {
							viewer.data().add_edges(w.row(i), w.row(i + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
						}
						if (wireinter && (!Intersection_al))
						{
							Is_Intersection[step / 3] = 1;
							Eigen::RowVector3d b2 = w.row(step / 3) - w.row((step / 3) + 1);
							double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
							double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
							double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
							double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
							double angle2 = acos(cos_angle2); // 夹角的弧度值
							double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
							if (angles[step / 3] > 0)
							{
								Intersection_angle[step / 3] = angle_degrees2;
							}
							else {
								Intersection_angle[step / 3] = -angle_degrees2;
							}
							Intersection_al = 1;
						}
						if (wireinter2 && (!cal_angle))
						{
							cal_angle = 1;
							V2_copy = V2; F2_copy = F2; V3_copy = V3; F3_copy = F3; w_copy = w; theta_copy = theta / 100;
							Eigen::RowVector3d pos2_copy = w.row(step / 3);
							w_copy(step / 3, 0) = w_copy((step / 3) + 1, 0) + (distances[step / 3] * (pos2_copy.x() - w_copy((step / 3) + 1, 0)) / sqrt(pow(pos2_copy.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos2_copy.y() - w_copy((step / 3) + 1, 1), 2)));
							w_copy(step / 3, 1) = w_copy((step / 3) + 1, 1) + (distances[step / 3] * (pos2_copy.y() - w_copy((step / 3) + 1, 1)) / sqrt(pow(pos2_copy.x() - w_copy((step / 3) + 1, 0), 2) + pow(pos2_copy.y() - w_copy((step / 3) + 1, 1), 2)));
						}
					}
					if (step / 3 > 0)
					{
						double angle_test = calculateAngle(w.row((step / 3) - 1), w.row(step / 3), w.row((step / 3) + 1));
						//std::cout << "angle_test: " << angle_test << std::endl;
					}

					//viewer.data().add_edges(w.row(step / 3), w.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
					// Check for collisions
					bool collision = igl::copyleft::cgal::intersect_other(V1, F1, V2, F2, false, IF);
					// If a collision is detected, reverse the rotation direction
					if (collision) {
						theta = -theta;
					}
					//
					b = w.row(step / 3) - w.row((step / 3) + 1);
					double dot_product = a.x() * b.x() + a.y() * b.y(); // 向量点乘
					double a_length = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
					double b_length = sqrt(b.x() * b.x() + b.y() * b.y()); // 向量 cd 的长度
					double cos_angle = dot_product / (a_length * b_length); // 夹角的余弦值
					double angle = acos(cos_angle); // 夹角的弧度值
					double angle_degrees = angle * 180 / M_PI; // 将弧度值转换为角度值
					if (angle_degrees >= 130)
					{
						Is_Intersection[step / 3] = 1;
						viewer.data().add_edges(w.row(step / 3), w.row((step / 3) + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
					}
					else
					{
						viewer.data().add_edges(w.row(step / 3), w.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
					}
					if (angle_degrees >= abs(angles[step / 3]))
					{
						//step++;
						intersection = 0;
						Intersection_al = 0;
						threebending++;
					}
				}
				else
				{
					if (cal_angle)
					{
						if (!cal_angle_init)
						{
							//viewer.data_list[1].show_faces = false;
							//viewer.data_list[2].show_faces = false;
							//viewer.data_list[3].show_faces = true;
							//viewer.data_list[4].show_faces = true;
							viewer.data_list[3].set_vertices(V2_copy);
							viewer.data_list[4].set_vertices(V3_copy);
							cal_angle_init = 1;
							//viewer.data().clear_edges();
							for (int i = 0; i < step / 3; i++)
							{
								bool wireinter = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V1, F1, hits2) && ((hits2[0].t < 1) && (hits2[0].t > 0));
								bool wireinter2 = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V2_copy, F2_copy, hits3) && ((hits3[0].t < 1) && (hits3[0].t > 0));
								if (wireinter || wireinter2)
								{
									//viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
								}
								else {
									//viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
								}
							}
							//viewer.data().add_edges(w_copy.row(step / 3), w_copy.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
						}
						if (cal_angle_init)
						{
							//viewer.data().clear_edges();
							pos_copy = w_copy.row((step / 3) - 1);
							w_copy((step / 3) - 1, 0) = (pos_copy.x() - w_copy((step / 3), 0)) * cos(theta_copy) + (pos_copy.y() - w_copy((step / 3), 1)) * sin(theta_copy) + w_copy((step / 3), 0);
							w_copy((step / 3) - 1, 1) = (pos_copy.y() - w_copy((step / 3), 1)) * cos(theta_copy) - (pos_copy.x() - w_copy((step / 3), 0)) * sin(theta_copy) + w_copy((step / 3), 1);
							for (int j = 0; j < ((step / 3) - 1); j++)
							{
								Eigen::RowVector3d p2_rotated = compute_rotation_point(pos, w_copy.row((step / 3) - 1), w_copy.row(j));
								w_copy.row(j) = p2_rotated;
							}
							pos2 = pos;
							for (int i = 0; i < step / 3; i++)
							{
								bool wireinter = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V1, F1, hits2) && ((hits2[0].t < 1) && (hits2[0].t > 0));
								bool wireinter2 = igl::ray_mesh_intersect(w_copy.row(i + 1), w_copy.row(i) - w_copy.row(i + 1), V2_copy, F2_copy, hits3) && ((hits3[0].t < 1) && (hits3[0].t > 0));
								if (wireinter || wireinter2)
								{
									//viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(1.0, 0.0, 0.0));
								}
								else {
									//viewer.data().add_edges(w_copy.row(i), w_copy.row(i + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
								}
								if (wireinter2)
								{
									cal_angle_ok = 1;
								}
							}
							//viewer.data().add_edges(w_copy.row(step / 3), w_copy.row((step / 3) + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
							//viewer.core().is_animating = !viewer.core().is_animating;
							if (cal_angle_ok == 0)
							{
								//计算角度
								Eigen::RowVector3d b2 = w_copy.row((step / 3)) - w_copy.row((step / 3) - 1);
								double dot_product2 = a.x() * b2.x() + a.y() * b2.y(); // 向量点乘
								double a_length2 = sqrt(a.x() * a.x() + a.y() * a.y()); // 向量 ab 的长度
								double b_length2 = sqrt(b2.x() * b2.x() + b2.y() * b2.y()); // 向量 cd 的长度
								double cos_angle2 = dot_product2 / (a_length2 * b_length2); // 夹角的余弦值
								double angle2 = acos(cos_angle2); // 夹角的弧度值
								double angle_degrees2 = angle2 * 180 / M_PI; // 将弧度值转换为角度值
								if (angles[(step) / 3] > 0)
								{
									std::cout << "chang angle to:" << 180 - angle_degrees2 << std::endl;
								}
								else {
									{std::cout << "chang angle to:" << -(180 - angle_degrees2) << std::endl; }
								}
								step++;
								threebending = 0;
								cal_angle_init = 0;
								cal_angle = 0;
								//viewer.data_list[1].show_faces = true;
								//viewer.data_list[2].show_faces = true;
								//viewer.data_list[3].show_faces = false;
								//viewer.data_list[4].show_faces = false;
								//viewer.data().clear_edges();
								//for (int j = 0; j <= step / 3; j++)
								//{
									//viewer.data().add_edges(w.row(j), w.row(j + 1), Eigen::RowVector3d(0.0, 1.0, 0.0));
								//}
							}
							cal_angle_ok = 0;//
						}
					}
					else {
						step++;
						threebending = 0;
					}
				}
			}
			drawpic();

		}
		return false;
	};

	// 定义用于启动和停止旋转的键盘回调功能
	viewer.callback_key_down = &key_down;
	// 启动渲染循环
	viewer.launch();
}


Collision_check::Collision_check() 
	: a(30, 0, 0), wireStartPoint(498.2, -240.48, 185.0), center3(492.2, -240.46, 168.8) {}

double Collision_check::angle_cal(const Eigen::RowVector3d& p1, const Eigen::RowVector3d& p2, const Eigen::RowVector3d& p3)
{
	// Calculate vectors between points
	Eigen::Vector3d vector1 = p1 - p2;
	Eigen::Vector3d vector2 = p3 - p2;
	// Calculate dot product of the two vectors
	double dotProduct = vector1.dot(vector2);
	// Calculate magnitudes of the vectors
	double magnitude1 = vector1.norm();
	double magnitude2 = vector2.norm();
	// Calculate the cosine of the angle between the vectors
	double cosine = dotProduct / (magnitude1 * magnitude2);
	// Calculate the angle in radians
	double angleInRadians = std::acos(cosine);
	// Convert angle from radians to degrees
	double angleInDegrees = angleInRadians * 180.0 / M_PI;

	return angleInDegrees;
}

bool Collision_check::isIntersecting(const Eigen::Vector3d& point_start, const Eigen::Vector3d& point_end, const Eigen::RowVector3d& min_corner, const Eigen::RowVector3d& max_corner)
{
	// 检查线段的起点和终点是否在包围盒内或与包围盒相交
	if (point_start[0] >= min_corner[0] && point_start[0] <= max_corner[0] &&
		point_start[1] >= min_corner[1] && point_start[1] <= max_corner[1] &&
		point_start[2] >= min_corner[2] && point_start[2] <= max_corner[2]) {
		return true; // 起点在包围盒内
	}

	if (point_end[0] >= min_corner[0] && point_end[0] <= max_corner[0] &&
		point_end[1] >= min_corner[1] && point_end[1] <= max_corner[1] &&
		point_end[2] >= min_corner[2] && point_end[2] <= max_corner[2]) {
		return true; // 终点在包围盒内
	}


	return false;


	////检测一个线段point_start-point_end是否与一个长方体相交
	////min_corner和max_corner分别是长方体的最小和最大顶点
 //   //如果相交，返回true，否则返回false
	//for (int i = 0; i < 3; i++) {
	//	if (point_start[i] > max_corner[i] && point_end[i] > max_corner[i]) {
	//		return false; // 线段在AABB的外部
	//	}
	//	if (point_start[i] < min_corner[i] && point_end[i] < min_corner[i]) {
	//		return false; // 线段在AABB的外部
	//	}
	//}
	//return true;
}

bool Collision_check::isIntersecting2(const Eigen::Vector3d& point_start, const Eigen::Vector3d& point_end, const Eigen::RowVector3d& min_corner_, const Eigen::RowVector3d& max_corner_)
{

	// 检查线段的起点和终点是否在包围盒内或与包围盒相交
	if (point_start[0] >= min_corner[0] && point_start[0] <= max_corner[0] &&
		point_start[1] >= min_corner[1] && point_start[1] <= max_corner[1] &&
		point_start[2] >= min_corner[2] && point_start[2] <= max_corner[2]) {
		return true; // 起点在包围盒内
	}

	if (point_end[0] >= min_corner[0] && point_end[0] <= max_corner[0] &&
		point_end[1] >= min_corner[1] && point_end[1] <= max_corner[1] &&
		point_end[2] >= min_corner[2] && point_end[2] <= max_corner[2]) {
		return true; // 终点在包围盒内
	}

	
	return false;

}

bool areSegmentsIntersecting(const Point_3& a, const Point_3& b, const Point_3& c, const Point_3& d) {
	Segment_3 segment_ab(a, b);
	Segment_3 segment_cd(c, d);

	// 使用CGAL提供的函数判断线段是否相交
	return CGAL::do_intersect(segment_ab, segment_cd);
}

//3：beam除cro碰撞检测时间    2：cro碰撞检测时间 其他碰撞检测时间
bool Collision_check::checkCollision(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx) {
	
	Eigen::MatrixXd w(anglesNumber, 3);
	for (int i = 0; i < w.rows(); ++i) {
		w.row(i) = wireStartPoint;
	}

	igl::readPLY("..\\file\\PLY\\cube.PLY", V1, F1);
	igl::readPLY("..\\file\\PLY\\rack_pin.PLY", V2, F2);
	igl::readPLY("..\\file\\PLY\\bender_gear.PLY", V3, F3);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\cube.PLY", V1, F1);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\rack_pin.PLY", V2, F2);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\bender_gear.PLY", V3, F3);

	min_corner = V1.colwise().minCoeff();
	max_corner = V1.colwise().maxCoeff();
	min_corner[0] -= 2; min_corner[1] -= 2; min_corner[2] -= 2;
	max_corner[0] += 2; max_corner[1] += 2; max_corner[2] += 2;
	step = 0;
	for (int bendingnum = startIdx; bendingnum < endIdx; bendingnum++)
	{
		//std::cout<<"bendingnum:"<<bendingnum<<std::endl;
		for (int i = 0; i <= step ; i++)
		{
			w(i, 0) += lengths[bendingnum];
			//std::cout<<"w("<<i<<", 0):"<<w(i, 0)<<"a" << std::endl;
		}


		if (abs(angles[bendingnum]) < 15|| (angles[bendingnum] > 114) || (angles[bendingnum] < -150))
		{
			return true;
		}
		/*if (abs(angles[bendingnum]) > 150 )
		{
			return true;
		}*/

		double angleInRadians;
		if (angles[bendingnum] > 0)
		{
			angleInRadians = -1.0 * M_PI / 380.0;
			//angleInRadians = -0.05;
		}
		else {
			angleInRadians = 1.0 * M_PI / 380.0;
			//angleInRadians = 0.05;
		}
		//for (int j = 0; j < abs(angles[bendingnum]); j++)
		//{
		//	
		//	Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
		//	for (int i = 0; i <= step; i++)
		//	{
		//		Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
		//		Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
		//		w.row(i) = rotatedPosition.transpose() + wireStartPoint;
		//		
		//	}
		//	for (int i = 0; i < step; i++)
		//	{
		//		//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
		//		bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
		//		if (wireinter)
		//		{
		//			return true;
		//		}
		//	}
		//}
		double angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
		while (1)
		{

			Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
			for (int i = 0; i <= step; i++)
			{
				Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
				Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
				w.row(i) = rotatedPosition.transpose() + wireStartPoint;

			}
			for (int i = 0; i < step; i++)
			{
				//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
				bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
				if (wireinter)
				{
					return true;
				}
			}
			angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
			if (angle_already >= abs(angles[bendingnum]))
			{
				break;
			}
		}

		/*for (int i = 0; i <= step; i++)
		{
			std::cout << "w(" << i << ", 0):" << w(i, 0) << std::endl;
		}*/
		step++;
	}
	return false;
}
bool Collision_check::checkCollision2(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx) {
	auto start_CC = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd w(anglesNumber, 3);
	for (int i = 0; i < w.rows(); ++i) {
		w.row(i) = wireStartPoint;
	}

	igl::readPLY("..\\file\\PLY\\cube.PLY", V1, F1);
	igl::readPLY("..\\file\\PLY\\rack_pin.PLY", V2, F2);
	igl::readPLY("..\\file\\PLY\\bender_gear.PLY", V3, F3);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\cube.PLY", V1, F1);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\rack_pin.PLY", V2, F2);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\bender_gear.PLY", V3, F3);

	min_corner = V1.colwise().minCoeff();
	max_corner = V1.colwise().maxCoeff();
	min_corner[0] -= 2; min_corner[1] -= 2; min_corner[2] -= 2;
	max_corner[0] += 2; max_corner[1] += 2; max_corner[2] += 2;
	step = 0;
	for (int bendingnum = startIdx; bendingnum < endIdx; bendingnum++)
	{
		//std::cout<<"bendingnum:"<<bendingnum<<std::endl;
		for (int i = 0; i <= step; i++)
		{
			w(i, 0) += lengths[bendingnum];
			//std::cout<<"w("<<i<<", 0):"<<w(i, 0)<<"a" << std::endl;
		}

		//if (abs(angles[bendingnum]) < 15 || (angles[bendingnum] > 114) || (angles[bendingnum] < -150))  abs(angles[bendingnum]) > 150
		if (abs(angles[bendingnum]) < 15 || (angles[bendingnum] > 114) || (angles[bendingnum] < -150))
		{
			auto end_CC = std::chrono::high_resolution_clock::now();
			timeCC2 += end_CC - start_CC;
			return true;
		}

		double angleInRadians;
		if (angles[bendingnum] > 0)
		{
			angleInRadians = -1.0 * M_PI / 380.0;
			//angleInRadians = -0.05;
		}
		else {
			angleInRadians = 1.0 * M_PI / 380.0;
			//angleInRadians = 0.05;
		}
		//for (int j = 0; j < abs(angles[bendingnum]); j++)
		//{
		//	
		//	Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
		//	for (int i = 0; i <= step; i++)
		//	{
		//		Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
		//		Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
		//		w.row(i) = rotatedPosition.transpose() + wireStartPoint;
		//		
		//	}
		//	for (int i = 0; i < step; i++)
		//	{
		//		//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
		//		bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
		//		if (wireinter)
		//		{
		//			return true;
		//		}
		//	}
		//}
		double angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
		while (1)
		{

			Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
			for (int i = 0; i <= step; i++)
			{
				Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
				Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
				w.row(i) = rotatedPosition.transpose() + wireStartPoint;

			}
			for (int i = 0; i < step; i++)
			{
				//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
				bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
				if (wireinter)
				{
					auto end_CC2 = std::chrono::high_resolution_clock::now();
					timeCC2 += end_CC2 - start_CC;
					return true;
				}
			}
			angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
			if (angle_already >= abs(angles[bendingnum]))
			{
				break;
			}
		}

		/*for (int i = 0; i <= step; i++)
		{
			std::cout << "w(" << i << ", 0):" << w(i, 0) << std::endl;
		}*/
		step++;
	}
	auto end_CC3 = std::chrono::high_resolution_clock::now();
	timeCC2 += end_CC3 - start_CC;
	return false;
}
bool Collision_check::checkCollision3(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx) {
	auto start_CC = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd w(anglesNumber, 3);
	for (int i = 0; i < w.rows(); ++i) {
		w.row(i) = wireStartPoint;
	}

	igl::readPLY("..\\file\\PLY\\cube.PLY", V1, F1);
	igl::readPLY("..\\file\\PLY\\rack_pin.PLY", V2, F2);
	igl::readPLY("..\\file\\PLY\\bender_gear.PLY", V3, F3);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\cube.PLY", V1, F1);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\rack_pin.PLY", V2, F2);
	//igl::readPLY("E:\\WireArtPjt\\model\\PLY\\bender_gear.PLY", V3, F3);

	min_corner = V1.colwise().minCoeff();
	max_corner = V1.colwise().maxCoeff();
	min_corner[0] -= 2; min_corner[1] -= 2; min_corner[2] -= 2;
	max_corner[0] += 2; max_corner[1] += 2; max_corner[2] += 2;

	step = 0;
	for (int bendingnum = startIdx; bendingnum < endIdx; bendingnum++)
	{
		//std::cout<<"bendingnum:"<<bendingnum<<std::endl;
		for (int i = 0; i <= step; i++)
		{
			w(i, 0) += lengths[bendingnum];
			//std::cout<<"w("<<i<<", 0):"<<w(i, 0)<<"a" << std::endl;
		}

		//if (abs(angles[bendingnum]) < 15 || (angles[bendingnum] > 114) || (angles[bendingnum] < -150))  abs(angles[bendingnum]) > 150
		if (abs(angles[bendingnum]) < 15 || (angles[bendingnum] > 114) || (angles[bendingnum] < -150))
		{
			auto end_CC = std::chrono::high_resolution_clock::now();
			timeCC1 += end_CC - start_CC;
			return true;
		}

		double angleInRadians;
		if (angles[bendingnum] > 0)
		{
			angleInRadians = -1.0 * M_PI / 380.0;
			//angleInRadians = -0.05;
		}
		else {
			angleInRadians = 1.0 * M_PI / 380.0;
			//angleInRadians = 0.05;
		}
		//for (int j = 0; j < abs(angles[bendingnum]); j++)
		//{
		//	
		//	Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
		//	for (int i = 0; i <= step; i++)
		//	{
		//		Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
		//		Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
		//		w.row(i) = rotatedPosition.transpose() + wireStartPoint;
		//		
		//	}
		//	for (int i = 0; i < step; i++)
		//	{
		//		//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
		//		bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
		//		if (wireinter)
		//		{
		//			return true;
		//		}
		//	}
		//}
		double angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
		while (1)
		{

			Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
			for (int i = 0; i <= step; i++)
			{
				Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
				Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
				w.row(i) = rotatedPosition.transpose() + wireStartPoint;

			}
			for (int i = 0; i < step; i++)
			{
				//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
				bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
				if (wireinter)
				{
					auto end_CC2 = std::chrono::high_resolution_clock::now();
					timeCC1 += end_CC2 - start_CC;
					return true;
				}
			}
			angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
			if (angle_already >= abs(angles[bendingnum]))
			{
				break;
			}
		}

		/*for (int i = 0; i <= step; i++)
		{
			std::cout << "w(" << i << ", 0):" << w(i, 0) << std::endl;
		}*/
		step++;
	}
	auto end_CC3 = std::chrono::high_resolution_clock::now();
	timeCC1 += end_CC3 - start_CC;
	return false;
}

Eigen::RowVector2d Collision_check::calculatePinCenter(const Eigen::RowVector2d& arcCenter, const Eigen::RowVector2d& wireStartPoint, double pinDistance, double radius) {
	// 计算 wireStartPoint 到 arcCenter 的向量
	Eigen::RowVector2d wireToArc = arcCenter - wireStartPoint;

	// 计算 wireToArc 的单位向量
	Eigen::RowVector2d wireToArcUnit = wireToArc.normalized();

	// 计算垂直于 wireToArc 的单位向量
	Eigen::RowVector2d perpendicularUnit(-wireToArcUnit(1), wireToArcUnit(0));

	// 计算 pinCenter 的位置
	Eigen::RowVector2d pinCenter = wireStartPoint + pinDistance * wireToArcUnit + radius * perpendicularUnit;

	return pinCenter;
}

Eigen::RowVector3d Collision_check::calculateCircleCenter(const Eigen::RowVector3d& startPoint, const Eigen::RowVector3d& finalPoint, double curvature_now) {
	
	//计算起始点到终止点的向量
	Eigen::RowVector3d directionVector = finalPoint - startPoint;

	//计算vector的垂直向量
	Eigen::RowVector3d perpendicularVector(-directionVector.y(), directionVector.x(), directionVector.z());
	Eigen::RowVector3d normalVector = perpendicularVector.normalized();

	// 计算圆弧的中点
	Eigen::RowVector3d midPoint = (startPoint + finalPoint) / 2.0;

	//计算finalPoint  startPoint的距离长度
	double dis = (finalPoint - startPoint).norm();
	double radius = abs(1.0 / curvature_now);
	double findis= sqrt(radius * radius - dis * dis / 4.0);
	
	Eigen::RowVector3d center;

	// 计算圆心坐标
	if (curvature_now > 0.0) {//顺时针
		center = midPoint - findis * normalVector;
	}
	else
	{
		center = midPoint + findis * normalVector;
	}
	
	return center;

}

//逆时针
Eigen::RowVector3d Collision_check::moveOnSphere(const Eigen::RowVector3d& point, const Eigen::RowVector3d& center, double arcLength, double curvature_now) {
	if (curvature_now == 0.0) {
		std::cout<<"curvature_now == 0.0"<<std::endl;
		return point;
	}
	double radius = 1 / abs(curvature_now);
	double startAngle = std::atan2(point(1) - center(1), point(0) - center(0));
    double deltaTheta;
	// 计算点在运动过程中绕圆心的旋转角度
	if(curvature_now>0)//顺时针
	{
		deltaTheta = -arcLength / radius;
	}
	else//逆时针
	{
		deltaTheta = arcLength / radius;
	}
	double finalAngle = startAngle + deltaTheta;
	// 将旋转后的位置向量重新加到运动圆心上，得到新的点坐标
	Eigen::RowVector3d newPoint;
	newPoint(0) = center(0) + radius * std::cos(finalAngle);
	newPoint(1) = center(1) + radius * std::sin(finalAngle);
	newPoint(2) = point(2);
	return newPoint;
}

bool Collision_check::checkCollision4(const std::vector<double>& lengths, const std::vector<double>& angles, int startIdx, int endIdx, const std::vector<bool>& segArc, const std::vector<double>& curvature) {

	Eigen::MatrixXd w(anglesNumber, 3);
	for (int i = 0; i < w.rows(); ++i) {
		w.row(i) = wireStartPoint;
	}

	igl::readPLY("..\\file\\PLY\\cube.PLY", V1, F1);
	igl::readPLY("..\\file\\PLY\\rack_pin.PLY", V2, F2);
	igl::readPLY("..\\file\\PLY\\bender_gear.PLY", V3, F3);

	min_corner = V1.colwise().minCoeff();
	max_corner = V1.colwise().maxCoeff();

	//Eigen::RowVector3d pinCenter(487.79, -213.64, 189.41); //wireStartPoint(498.2, -240.48, 185.0)
	double pinRadius = 6.2;//半径
	double pinDistance = 28.8;//pin距离出线口距离

	step = 0;
	for (int bendingnum = startIdx; bendingnum < endIdx; bendingnum++)
	{
		if(segArc[bendingnum] == false)
		{
			//std::cout<<"bendingnum:"<<bendingnum<<std::endl;
			for (int i = 0; i <= step; i++)
			{
				w(i, 0) += lengths[bendingnum];
				//std::cout<<"w("<<i<<", 0):"<<w(i, 0)<<"a" << std::endl;
			}

			if (abs(angles[bendingnum]) < 15 || (angles[bendingnum] > 114) || (angles[bendingnum] < -150))
			{
				return true;
			}
			//if (abs(angles[bendingnum]) > 150)
			//{
			//	return true;
			//}

			double angleInRadians;
			if (angles[bendingnum] > 0)
			{
				angleInRadians = -1.0 * M_PI / 380.0;
			}
			else {
				angleInRadians = 1.0 * M_PI / 380.0;
			}
			double angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
			while (1)
			{

				Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
				for (int i = 0; i <= step; i++)
				{
					Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
					Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
					w.row(i) = rotatedPosition.transpose() + wireStartPoint;

				}
				for (int i = 0; i < step; i++)
				{
					////bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
					//bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
					//if (wireinter)
					//{
					//	return true;
					//}
					Eigen::RowVector3d startPoint = w.row(i + 1);
					Eigen::RowVector3d finalPoint = w.row(i);
					if (segArc[startIdx + i])
					{
						Eigen::RowVector3d firstPoint = w.row(i+1);
						Eigen::RowVector3d secondPoint = w.row(i+1);
						double curvature_now = curvature[startIdx + i];//曲率
						double radius_now = 1 / abs(curvature_now);
						double arc_length_now = lengths[startIdx + i]; //弧长
						double angle_now = arc_length_now / radius_now; //弧度
						double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

						Eigen::RowVector3d arcCenter_now = calculateCircleCenter(startPoint, finalPoint, curvature_now);

						for (double k = 0; k <arc_length_now; k++)
						{
							firstPoint = moveOnSphere(firstPoint, arcCenter_now, 1, curvature_now);
							//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
							bool wireinter = isIntersecting(firstPoint, secondPoint, min_corner, max_corner);
							if (wireinter)
							{
								return true;
							}
							secondPoint = firstPoint;
						}
					}
					else
					{
						//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
						bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
						if (wireinter)
						{
							return true;
						}
					}
				}
				/*for (int m = 0; m <= step; m++)
				{
					Eigen::RowVector3d segarcstart1 = w.row(m+1);
					Eigen::RowVector3d segarcend1 = w.row(m);
					for (int n = m + 1; n <= step; n++)
					{
						Eigen::RowVector3d segarcstart2 = w.row(n + 1);
						Eigen::RowVector3d segarcend2 = w.row(n);
						Point_3 a(segarcstart1(0), segarcstart1(1), segarcstart1(2));
						Point_3 b(segarcend1(0), segarcend1(1), segarcend1(2));
						Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
						Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
						if (areSegmentsIntersecting(a, b, c, d)) {
							return true;
						}				
					}
				}*/
				/*for (int m = 0; m <= step; m++)
				{
					Eigen::RowVector3d segarcstart1 = w.row(m + 1);
					Eigen::RowVector3d segarcend1 = w.row(m);
					for (int n = m + 1; n <= step; n++)
					{
						Eigen::RowVector3d segarcstart2 = w.row(n + 1);
						Eigen::RowVector3d segarcend2 = w.row(n);
						if ((!segArc[startIdx + m]) && (!segArc[startIdx + n]))//false线段
						{

							Point_3 a(segarcstart1(0), segarcstart1(1), segarcstart1(2));
							Point_3 b(segarcend1(0), segarcend1(1), segarcend1(2));
							Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
							Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
							if (areSegmentsIntersecting(a, b, c, d)) {
								return true;
							}
						}
						else if ((segArc[startIdx + m]) && (segArc[startIdx + n]))//圆弧
						{
							Eigen::RowVector3d firstPoint1 = w.row(m + 1);
							Eigen::RowVector3d secondPoint1 = w.row(m + 1);
							double curvature_now1 = curvature[startIdx + m];//曲率
							double radius_now1 = 1 / abs(curvature_now1);
							double arc_length_now1 = lengths[startIdx + m]; //弧长
							double angle_now1 = arc_length_now1 / radius_now1; //弧度
							double angle_now_degrees1 = angle_now1 * (180.0 / M_PI); //角度
							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart1, segarcend1, curvature_now1);

							Eigen::RowVector3d firstPoint2 = w.row(n + 1);
							Eigen::RowVector3d secondPoint2 = w.row(n + 1);
							double curvature_now2 = curvature[startIdx + n];//曲率
							double radius_now2 = 1 / abs(curvature_now2);
							double arc_length_now2 = lengths[startIdx + n]; //弧长
							double angle_now2 = arc_length_now2 / radius_now2; //弧度
							double angle_now_degrees2 = angle_now2 * (180.0 / M_PI); //角度
							Eigen::RowVector3d arcCenter_now2 = calculateCircleCenter(segarcstart2, segarcend2, curvature_now2);
							for (double k1 = 0; k1 < angle_now_degrees1; k1++)
							{
								for (double k2 = 0; k2 < angle_now_degrees2; k2++)
								{
									firstPoint1 = moveOnSphere(firstPoint1, arcCenter_now1, radius_now1 * M_PI / 180.0, curvature_now1);
									firstPoint2 = moveOnSphere(firstPoint2, arcCenter_now2, radius_now2 * M_PI / 180.0, curvature_now2);
									Point_3 a(firstPoint1(0), firstPoint1(1), firstPoint1(2));
									Point_3 b(secondPoint1(0), secondPoint1(1), secondPoint1(2));
									Point_3 c(firstPoint2(0), firstPoint2(1), firstPoint2(2));
									Point_3 d(secondPoint2(0), secondPoint2(1), secondPoint2(2));
									if (areSegmentsIntersecting(a, b, c, d)) {
										return true;
									}
									secondPoint1 = firstPoint1;
									secondPoint2 = firstPoint2;
								}
							}
						}
						else if ((segArc[startIdx + m]) && (!segArc[startIdx + n]))//圆弧 线段
						{
							Eigen::RowVector3d firstPoint = w.row(m + 1);
							Eigen::RowVector3d secondPoint = w.row(m + 1);
							double curvature_now = curvature[startIdx + m];//曲率
							double radius_now = 1 / abs(curvature_now);
							double arc_length_now = lengths[startIdx + m]; //弧长
							double angle_now = arc_length_now / radius_now; //弧度
							double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart1, segarcend1, curvature_now);

							for (double k = 0; k < angle_now_degrees; k++)
							{
								firstPoint = moveOnSphere(firstPoint, arcCenter_now1, radius_now * M_PI / 180.0, curvature_now);
								//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
								Point_3 a(firstPoint(0), firstPoint(1), firstPoint(2));
								Point_3 b(secondPoint(0), secondPoint(1), secondPoint(2));
								Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
								Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
								if (areSegmentsIntersecting(a, b, c, d)) {
									return true;
								}
								secondPoint = firstPoint;
							}
						}
						else if ((!segArc[startIdx + m]) && (segArc[startIdx + n]))
						{
							Eigen::RowVector3d firstPoint = w.row(n + 1);
							Eigen::RowVector3d secondPoint = w.row(n + 1);
							double curvature_now = curvature[startIdx + n];//曲率
							double radius_now = 1 / abs(curvature_now);
							double arc_length_now = lengths[startIdx + n]; //弧长
							double angle_now = arc_length_now / radius_now; //弧度
							double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart2, segarcend2, curvature_now);

							for (double k = 0; k < angle_now_degrees; k++)
							{
								firstPoint = moveOnSphere(firstPoint, arcCenter_now1, radius_now * M_PI / 180.0, curvature_now);
								//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
								Point_3 a(firstPoint(0), firstPoint(1), firstPoint(2));
								Point_3 b(secondPoint(0), secondPoint(1), secondPoint(2));
								Point_3 c(segarcstart1(0), segarcstart1(1), segarcstart1(2));
								Point_3 d(segarcend1(0), segarcend1(1), segarcend1(2));
								if (areSegmentsIntersecting(a, b, c, d)) {
									return true;
								}
								secondPoint = firstPoint;
							}
						}
					}
				}*/
				angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
				if (angle_already >= abs(angles[bendingnum]))
				{
					break;
				}
			}
		}
		else {
			if (curvature[bendingnum] == 0)
			{
				std::cout<<"bug: arc curvature = 0"<<std::endl;
				return false;
			}
			double radius = 1 / abs(curvature[bendingnum]);
			Eigen::RowVector3d arcCenter = wireStartPoint;
			if(curvature[bendingnum]>0)
			{arcCenter(1)=arcCenter(1)-radius;}
			else
			{arcCenter(1)=arcCenter(1)+radius;}
			//Eigen::RowVector3d pinCenter;
			//pinCenter(2) = wireStartPoint(2);
			//曲率为正，圆心在下方  另一方向不一定对
			//Eigen::RowVector2d arcCenter2d(arcCenter(0), arcCenter(1));
			//Eigen::RowVector2d wireStartPoint2d(wireStartPoint(0), wireStartPoint(1));
			//Eigen::RowVector2d pinCenter2d;
			//pinCenter2d = calculatePinCenter(arcCenter2d, wireStartPoint2d, pinDistance, radius+ pinRadius);
			//pinCenter(0) = pinCenter2d(0);
			//pinCenter(1) = pinCenter2d(1);
			//Eigen::RowVector3d pinCounter;
			//Eigen::RowVector3d arcToPinUnit = (pinCenter - arcCenter).normalized();
			//pinCounter = pinCenter + pinRadius * arcToPinUnit;
			//while (!w.row(step).isApprox(pinCounter))
			//{
			//	std::cout<<"pinCounter:"<< w.row(step) <<std::endl;
			//	w.row(step) = moveOnSphere(wireStartPoint, arcCenter, 5);
			//	//一起旋转 是否对？
			//	for (int i = 0; i < step; i++)
			//	{
			//		w.row(i) = moveOnSphere(wireStartPoint, arcCenter, 5);
			//	}
			//	//碰撞检测
			//	for (int i = 0; i < step; i++)
			//	{
			//		//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
			//		bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
			//		if (wireinter)
			//		{
			//			std::cout << "D";
			//			return true;
			//		}
			//	}
			//}
			
			for(double forwaedLengthAleady = 1; forwaedLengthAleady <= lengths[bendingnum]; forwaedLengthAleady = forwaedLengthAleady + 1)
			{
				//一起旋转 是否对？
				for (int i = 0; i <= step; i++)
				{
					w.row(i) = moveOnSphere(w.row(i), arcCenter, 1, curvature[bendingnum]);
				}
				//碰撞检测
				for (int i = 0; i < step; i++)
				{
					Eigen::RowVector3d startPoint = w.row(i+1);
					Eigen::RowVector3d finalPoint = w.row(i);
					if (segArc[startIdx + i])
					{
						Eigen::RowVector3d firstPoint = w.row(i + 1);
						Eigen::RowVector3d secondPoint = w.row(i + 1);
						double curvature_now = curvature[startIdx + i];//曲率
						double radius_now = 1 / abs(curvature_now);
						double arc_length_now = lengths[startIdx + i]; //弧长
						double angle_now = arc_length_now / radius_now; //弧度
						double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

						Eigen::RowVector3d arcCenter_now = calculateCircleCenter(startPoint, finalPoint, curvature_now);

						for (double k = 0; k < arc_length_now; k++)
						{
							firstPoint = moveOnSphere(firstPoint, arcCenter_now, 1, curvature_now);
							//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
							bool wireinter = isIntersecting(firstPoint, secondPoint, min_corner, max_corner);
							if (wireinter)
							{
								return true;
							}
							secondPoint = firstPoint;
						}
						//curvature = segment["curvature"]
						/*arc_length = segment["arc_length"]
							//radius = 1 / abs(curvature)
							angle_change = math.degrees(arc_length / radius)  # 弧长转换为角度

							angle_rad = math.radians(abs(segment["angle"]))
							for arc_segment in range(int(angle_change)) :


								angle_rad += math.radians(1)  # 增加1度
								x += math.cos(angle_rad) * arc_length / angle_change
								if curvature < 0:
						y += math.sin(angle_rad) * arc_length / angle_change
								else:
						y -= math.sin(angle_rad) * arc_length / angle_change*/
						//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
						
					}
					else
					{
						//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
						bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
						if (wireinter)
						{
							return true;
						}
					}
					
				}
				/*for (int m = 0; m <= step; m++)
				{
					Eigen::RowVector3d segarcstart1 = w.row(m + 1);
					Eigen::RowVector3d segarcend1 = w.row(m);
					for (int n = m + 1; n <= step; n++)
					{
						Eigen::RowVector3d segarcstart2 = w.row(n + 1);
						Eigen::RowVector3d segarcend2 = w.row(n);
						Point_3 a(segarcstart1(0), segarcstart1(1), segarcstart1(2));
						Point_3 b(segarcend1(0), segarcend1(1), segarcend1(2));
						Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
						Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
						if (areSegmentsIntersecting(a, b, c, d)) {
							return true;
						}

					}
				}*/
				/*for (int m = 0; m <= step; m++)
				{
					Eigen::RowVector3d segarcstart1 = w.row(m + 1);
					Eigen::RowVector3d segarcend1 = w.row(m);
					for (int n = m + 1; n <= step; n++)
					{
						Eigen::RowVector3d segarcstart2 = w.row(n + 1);
						Eigen::RowVector3d segarcend2 = w.row(n);
						if ((!segArc[startIdx + m]) && (!segArc[startIdx + n]))//false线段
						{

							Point_3 a(segarcstart1(0), segarcstart1(1), segarcstart1(2));
							Point_3 b(segarcend1(0), segarcend1(1), segarcend1(2));
							Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
							Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
							if (areSegmentsIntersecting(a, b, c, d)) {
								return true;
							}
						}
						else if ((segArc[startIdx + m]) && (segArc[startIdx + n]))//圆弧
						{
							Eigen::RowVector3d firstPoint1 = w.row(m + 1);
							Eigen::RowVector3d secondPoint1 = w.row(m + 1);
							double curvature_now1 = curvature[startIdx + m];//曲率
							double radius_now1 = 1 / abs(curvature_now1);
							double arc_length_now1 = lengths[startIdx + m]; //弧长
							double angle_now1 = arc_length_now1 / radius_now1; //弧度
							double angle_now_degrees1 = angle_now1 * (180.0 / M_PI); //角度
							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart1, segarcend1, curvature_now1);

							Eigen::RowVector3d firstPoint2 = w.row(n + 1);
							Eigen::RowVector3d secondPoint2 = w.row(n + 1);
							double curvature_now2 = curvature[startIdx + n];//曲率
							double radius_now2 = 1 / abs(curvature_now2);
							double arc_length_now2 = lengths[startIdx + n]; //弧长
							double angle_now2 = arc_length_now2 / radius_now2; //弧度
							double angle_now_degrees2 = angle_now2 * (180.0 / M_PI); //角度
							Eigen::RowVector3d arcCenter_now2 = calculateCircleCenter(segarcstart2, segarcend2, curvature_now2);
							for (double k1 = 0; k1 < angle_now_degrees1; k1++)
							{
								for (double k2 = 0; k2 < angle_now_degrees2; k2++)
								{
									firstPoint1 = moveOnSphere(firstPoint1, arcCenter_now1, radius_now1 * M_PI / 180.0, curvature_now1);
									firstPoint2 = moveOnSphere(firstPoint2, arcCenter_now2, radius_now2 * M_PI / 180.0, curvature_now2);
									Point_3 a(firstPoint1(0), firstPoint1(1), firstPoint1(2));
									Point_3 b(secondPoint1(0), secondPoint1(1), secondPoint1(2));
									Point_3 c(firstPoint2(0), firstPoint2(1), firstPoint2(2));
									Point_3 d(secondPoint2(0), secondPoint2(1), secondPoint2(2));
									if (areSegmentsIntersecting(a, b, c, d)) {
										return true;
									}
									secondPoint1 = firstPoint1;
									secondPoint2 = firstPoint2;
								}
							}
						}
						else if ((segArc[startIdx + m]) && (!segArc[startIdx + n]))//圆弧 线段
						{
							Eigen::RowVector3d firstPoint = w.row(m + 1);
							Eigen::RowVector3d secondPoint = w.row(m + 1);
							double curvature_now = curvature[startIdx + m];//曲率
							double radius_now = 1 / abs(curvature_now);
							double arc_length_now = lengths[startIdx + m]; //弧长
							double angle_now = arc_length_now / radius_now; //弧度
							double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart1, segarcend1, curvature_now);

							for (double k = 0; k < angle_now_degrees; k++)
							{
								firstPoint = moveOnSphere(firstPoint, arcCenter_now1, radius_now * M_PI / 180.0, curvature_now);
								//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
								Point_3 a(firstPoint(0), firstPoint(1), firstPoint(2));
								Point_3 b(secondPoint(0), secondPoint(1), secondPoint(2));
								Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
								Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
								if (areSegmentsIntersecting(a, b, c, d)) {
									return true;
								}
								secondPoint = firstPoint;
							}
						}
						else if ((!segArc[startIdx + m]) && (segArc[startIdx + n]))
						{
							Eigen::RowVector3d firstPoint = w.row(n + 1);
							Eigen::RowVector3d secondPoint = w.row(n + 1);
							double curvature_now = curvature[startIdx + n];//曲率
							double radius_now = 1 / abs(curvature_now);
							double arc_length_now = lengths[startIdx + n]; //弧长
							double angle_now = arc_length_now / radius_now; //弧度
							double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart2, segarcend2, curvature_now);

							for (double k = 0; k < angle_now_degrees; k++)
							{
								firstPoint = moveOnSphere(firstPoint, arcCenter_now1, radius_now * M_PI / 180.0, curvature_now);
								//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
								Point_3 a(firstPoint(0), firstPoint(1), firstPoint(2));
								Point_3 b(secondPoint(0), secondPoint(1), secondPoint(2));
								Point_3 c(segarcstart1(0), segarcstart1(1), segarcstart1(2));
								Point_3 d(segarcend1(0), segarcend1(1), segarcend1(2));
								if (areSegmentsIntersecting(a, b, c, d)) {
									return true;
								}
								secondPoint = firstPoint;
							}
						}
					}
				}*/
                
				
			}
			if (abs(angles[bendingnum]) < 15 || (angles[bendingnum] > 114) || (angles[bendingnum] < -150))
			{
				return true;
			}
			//if (abs(angles[bendingnum]) > 150)
			//{
			//	return true;
			//}

			double angleInRadians;
			if (angles[bendingnum] > 0)
			{
				angleInRadians = -1.0 * M_PI / 380.0;
				//angleInRadians = -0.05;
			}
			else {
				angleInRadians = 1.0 * M_PI / 380.0;
				//angleInRadians = 0.05;
			}

			double angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
			while (1)
			{

				Eigen::Matrix3d rotationMatrix = Eigen::AngleAxisd(angleInRadians, Eigen::Vector3d::UnitZ()).matrix();
				for (int i = 0; i <= step; i++)
				{
					Eigen::Vector3d relativePosition = w.row(i).transpose() - wireStartPoint.transpose();
					Eigen::Vector3d rotatedPosition = rotationMatrix * relativePosition;
					w.row(i) = rotatedPosition.transpose() + wireStartPoint;

				}
				for (int i = 0; i < step; i++)
				{
					////bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
					//bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
					//if (wireinter)
					//{
					//	return true;
					//}
					Eigen::RowVector3d startPoint = w.row(i+ 1);
					Eigen::RowVector3d finalPoint = w.row(i );
					if (segArc[startIdx + i])
					{
						Eigen::RowVector3d firstPoint = w.row(i+1);
						Eigen::RowVector3d secondPoint = w.row(i+1);
						double curvature_now = curvature[startIdx + i];//曲率
						double radius_now = 1 / abs(curvature_now);
						double arc_length_now = lengths[startIdx + i]; //弧长
						double angle_now = arc_length_now / radius_now; //弧度
						double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

						Eigen::RowVector3d arcCenter_now = calculateCircleCenter(startPoint, finalPoint, curvature_now);

						for (double k = 0; k < arc_length_now; k++)
						{
							firstPoint = moveOnSphere(firstPoint, arcCenter_now, 1, curvature_now);
							//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
							bool wireinter = isIntersecting(firstPoint, secondPoint, min_corner, max_corner);
							if (wireinter)
							{
								return true;
							}
							secondPoint = firstPoint;
						}
						
					}
					else
					{
						//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
						bool wireinter = isIntersecting(w.row(i + 1), w.row(i), min_corner, max_corner);
						if (wireinter)
						{
							return true;
						}
					}

				}
				/*for (int m = 0; m <= step; m++)
				{
					Eigen::RowVector3d segarcstart1 = w.row(m + 1);
					Eigen::RowVector3d segarcend1 = w.row(m);
					for (int n = m + 1; n <= step; n++)
					{
						Eigen::RowVector3d segarcstart2 = w.row(n + 1);
						Eigen::RowVector3d segarcend2 = w.row(n);
						Point_3 a(segarcstart1(0), segarcstart1(1), segarcstart1(2));
						Point_3 b(segarcend1(0), segarcend1(1), segarcend1(2));
						Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
						Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
						if (areSegmentsIntersecting(a, b, c, d)) {
							return true;
						}

					}
				}*/
				/*for (int m = 0; m <= step; m++)
				{
					Eigen::RowVector3d segarcstart1 = w.row(m + 1);
					Eigen::RowVector3d segarcend1 = w.row(m);
					for (int n = m + 1; n <= step; n++)
					{
						Eigen::RowVector3d segarcstart2 = w.row(n + 1);
						Eigen::RowVector3d segarcend2 = w.row(n);
						if ((!segArc[startIdx + m]) && (!segArc[startIdx + n]))//false线段
						{

							Point_3 a(segarcstart1(0), segarcstart1(1), segarcstart1(2));
							Point_3 b(segarcend1(0), segarcend1(1), segarcend1(2));
							Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
							Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
							if (areSegmentsIntersecting(a, b, c, d)) {
								return true;
							}
						}
						else if ((segArc[startIdx + m]) && (segArc[startIdx + n]))//圆弧
						{
							Eigen::RowVector3d firstPoint1 = w.row(m + 1);
							Eigen::RowVector3d secondPoint1 = w.row(m + 1);
							double curvature_now1 = curvature[startIdx + m];//曲率
							double radius_now1 = 1 / abs(curvature_now1);
							double arc_length_now1 = lengths[startIdx + m]; //弧长
							double angle_now1 = arc_length_now1 / radius_now1; //弧度
							double angle_now_degrees1 = angle_now1 * (180.0 / M_PI); //角度
							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart1, segarcend1, curvature_now1);

							Eigen::RowVector3d firstPoint2 = w.row(n + 1);
							Eigen::RowVector3d secondPoint2 = w.row(n + 1);
							double curvature_now2 = curvature[startIdx + n];//曲率
							double radius_now2 = 1 / abs(curvature_now2);
							double arc_length_now2 = lengths[startIdx + n]; //弧长
							double angle_now2 = arc_length_now2 / radius_now2; //弧度
							double angle_now_degrees2 = angle_now2 * (180.0 / M_PI); //角度
							Eigen::RowVector3d arcCenter_now2 = calculateCircleCenter(segarcstart2, segarcend2, curvature_now2);
							for (double k1 = 0; k1 < angle_now_degrees1; k1++)
							{
								for (double k2 = 0; k2 < angle_now_degrees2; k2++)
								{
									firstPoint1 = moveOnSphere(firstPoint1, arcCenter_now1, radius_now1 * M_PI / 180.0, curvature_now1);
									firstPoint2 = moveOnSphere(firstPoint2, arcCenter_now2, radius_now2 * M_PI / 180.0, curvature_now2);
									Point_3 a(firstPoint1(0), firstPoint1(1), firstPoint1(2));
									Point_3 b(secondPoint1(0), secondPoint1(1), secondPoint1(2));
									Point_3 c(firstPoint2(0), firstPoint2(1), firstPoint2(2));
									Point_3 d(secondPoint2(0), secondPoint2(1), secondPoint2(2));
									if (areSegmentsIntersecting(a, b, c, d)) {
										return true;
									}
									secondPoint1 = firstPoint1;
									secondPoint2 = firstPoint2;
								}
							}
						}
						else if ((segArc[startIdx + m]) && (!segArc[startIdx + n]))//圆弧 线段
						{
							Eigen::RowVector3d firstPoint = w.row(m + 1);
							Eigen::RowVector3d secondPoint = w.row(m + 1);
							double curvature_now = curvature[startIdx + m];//曲率
							double radius_now = 1 / abs(curvature_now);
							double arc_length_now = lengths[startIdx + m]; //弧长
							double angle_now = arc_length_now / radius_now; //弧度
							double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart1, segarcend1, curvature_now);

							for (double k = 0; k < angle_now_degrees; k++)
							{
								firstPoint = moveOnSphere(firstPoint, arcCenter_now1, radius_now * M_PI / 180.0, curvature_now);
								//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
								Point_3 a(firstPoint(0), firstPoint(1), firstPoint(2));
								Point_3 b(secondPoint(0), secondPoint(1), secondPoint(2));
								Point_3 c(segarcstart2(0), segarcstart2(1), segarcstart2(2));
								Point_3 d(segarcend2(0), segarcend2(1), segarcend2(2));
								if (areSegmentsIntersecting(a, b, c, d)) {
									return true;
								}
								secondPoint = firstPoint;
							}
						}
						else if ((!segArc[startIdx + m]) && (segArc[startIdx + n]))
						{
							Eigen::RowVector3d firstPoint = w.row(n + 1);
							Eigen::RowVector3d secondPoint = w.row(n + 1);
							double curvature_now = curvature[startIdx + n];//曲率
							double radius_now = 1 / abs(curvature_now);
							double arc_length_now = lengths[startIdx + n]; //弧长
							double angle_now = arc_length_now / radius_now; //弧度
							double angle_now_degrees = angle_now * (180.0 / M_PI); //角度

							Eigen::RowVector3d arcCenter_now1 = calculateCircleCenter(segarcstart2, segarcend2, curvature_now);

							for (double k = 0; k < angle_now_degrees; k++)
							{
								firstPoint = moveOnSphere(firstPoint, arcCenter_now1, radius_now * M_PI / 180.0, curvature_now);
								//bool wireinter = igl::ray_mesh_intersect(w.row(i + 1), w.row(i) - w.row(i + 1), V1, F1, hits) && ((hits[0].t < 1) && (hits[0].t > 0));
								Point_3 a(firstPoint(0), firstPoint(1), firstPoint(2));
								Point_3 b(secondPoint(0), secondPoint(1), secondPoint(2));
								Point_3 c(segarcstart1(0), segarcstart1(1), segarcstart1(2));
								Point_3 d(segarcend1(0), segarcend1(1), segarcend1(2));
								if (areSegmentsIntersecting(a, b, c, d)) {
									return true;
								}
								secondPoint = firstPoint;
							}
						}
					}
				}*/
				angle_already = angle_cal(w.row(step), w.row(step + 1), w.row(step + 1) + a);
				if (angle_already >= abs(angles[bendingnum]))
				{
					break;
				}
			}
		}
		step++;
	}
	return false;
}

std::pair<int, int> findLongestNonCollidingSegment(const std::vector<double>& lengths, const std::vector<double>& angles, int startIndex) {
	
	int n = lengths.size();
	int startIdx = startIndex;
	int endIdx = startIndex;
	
	// Forward search
	for (int i = startIndex; i < (n - 1); ++i) {
		Collision_check collisionChecker;
		if (collisionChecker.checkCollision(lengths, angles, startIdx, endIdx + 1)) {
			break; // Collision occurred, stop forward search
		}
		endIdx++;
	}

	// Backward search
	for (int i = startIndex; i > 0; --i) {
		Collision_check collisionChecker;
		if (collisionChecker.checkCollision(lengths, angles, startIdx-1, endIdx)) {
			break; // Collision occurred, stop forward search
		}
		startIdx--;
	}
	return std::make_pair(startIdx, endIdx);
}

std::pair<int, int> findLongestNonCollidingSegment2(const std::vector<double>& lengths, const std::vector<double>& angles, int startIndex) {

	int n = lengths.size();
	int startIdx = startIndex;
	int endIdx = startIndex;
	
	// Backward search
	for (int i = startIndex; i > 0; --i) {
		Collision_check collisionChecker;
		if (collisionChecker.checkCollision(lengths, angles, startIdx - 1, endIdx)) {
			break; // Collision occurred, stop forward search
		}
		startIdx--;
	}

	return std::make_pair(startIdx, endIdx);
}