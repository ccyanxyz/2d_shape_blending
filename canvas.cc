#include <iostream>
#include <map>
#include <algorithm>
#include <Eigen/Dense>
#include <chrono>
#include <thread>
#include <QWidget>
#include <QTimer>
#include <QPainter>
#include <QBrush>
#include <QPixmap>
#include <QMouseEvent>
#include <QDebug>
#include <QMessageBox>
#include <QLabel>
#include <QPropertyAnimation>
#include <iostream>
#include <cmath>
#include <cassert>
#include <vector>
#include "canvas.h"

using namespace std;

void Canvas::clear()
{
	poly1.clear();
	poly2.clear();
	num = 0;
	for(int i = 0; i < 100; ++i) {
		for(int j = 0; j < 100; ++j) {
			sim[i][j] = 0;
		}
	}
	min_dist = 10000;
	mapping.clear();
	show_mapping = false;

	anchor_points.clear();
	show_anchor_points = false;

	trans_frames.clear();
	play_status = false;
}

void Canvas::draw_done()
{
	if(num <= 1) {
		num += 1;
		if(num == 2) {
			calc_poly_sim();
			calc_points_map();
			calc_best_affine_trans();
		}
		update();
	}
}

void Canvas::play()
{
	if(play_status == false) {
		trans_frames = interpolation();
		play_status = true;
	} else {
		play_status = false;
	}

	update();
}

void Canvas::paintEvent(QPaintEvent *)
{
	QPainter p(this);

	QBrush brush(Qt::black, Qt::SolidPattern);
	p.setPen(QColor(Qt::red));
	p.setBrush(brush);
	// paint control points
	for(size_t i = 0; i < poly1.size() && play_status == false; ++i){
		p.drawEllipse(poly1[i].x, poly1[i].y, 5, 5);
	}

	// paint line
	for(size_t i = 0; i < poly1.size() && play_status == false; ++i){
		if(i == poly1.size() - 1 && num >= 1) {
			p.drawLine(poly1[i].x, poly1[i].y, poly1[0].x, poly1[0].y);
		}

		if(i != poly1.size() - 1) {
			p.drawLine(poly1[i].x, poly1[i].y, poly1[i + 1].x, poly1[i + 1].y);
		}
	}

	for(size_t i = 0; i < poly2.size(); ++i) {
		p.drawEllipse(poly2[i].x, poly2[i].y, 5, 5);
	}

	for(size_t i = 0; i < poly2.size(); ++i) {
		if(i == poly2.size() - 1 && num >= 2) {
			p.drawLine(poly2[i].x, poly2[i].y, poly2[0].x, poly2[0].y);
		}

		if(i != poly2.size() - 1) {
			p.drawLine(poly2[i].x, poly2[i].y, poly2[i + 1].x, poly2[i + 1].y);
		}
	}

	// show points mapping
	if(show_mapping && mapping.size() == poly1.size() && play_status == false) {
		p.setPen(QColor(Qt::green));
		for(auto it = mapping.begin(); it != mapping.end(); ++it) {
			p.drawLine(poly1[it->first].x, poly1[it->first].y,\
					poly2[it->second].x, poly2[it->second].y);
		}
	}

	// show anchor points
	if(show_anchor_points && anchor_points.size() == 3 && play_status == false) {
		p.setPen(QColor(Qt::blue));
		for(auto i : anchor_points) {
			p.drawEllipse(poly1[i].x, poly1[i].y, 5, 5);
			p.drawEllipse(poly2[mapping[i]].x, poly2[mapping[i]].y, 5, 5);
			p.drawLine(poly1[i].x, poly1[i].y, poly2[mapping[i]].x,\
					poly2[mapping[i]].y);
		}
	}

	if(play_status && trans_frames.size() > 0) {
		if(pos >= trans_frames.size()) {
			play_status = false;
			pos = 0;
		} else {
			vector<Point> frame = trans_frames[pos];
			for(size_t i = 0; i < frame.size(); ++i) {
				p.drawLine(frame[i].x, frame[i].y, frame[(i + 1) % frame.size()].x,\
						frame[(i + 1) % frame.size()].y);
			}
			pos += 1;
			this_thread::sleep_for(chrono::milliseconds(25));
			update();
		}
	}
}

void Canvas::mousePressEvent(QMouseEvent *e)
{
	float x = e->pos().x();
	float y = e->pos().y();

	Point p(x, y);
	if(num == 0) {
		poly1.push_back(p);
	} else if(num == 1) {
		poly2.push_back(p);
	}
	cout << "x:" << x << ", y:" << y << endl;
	update();
}

double Canvas::calc_edge_len(Point &a, Point &b)
{
	double len = (a.x - b.x) * (a.x - b.x)\
				 + (a.y - b.y) * (a.y - b.y);

	return sqrt(len);
}

double Canvas::calc_angle(double a, double b, double c)
{
	double cosa = (b * b + c * c - a * a) / (2 * b * c);

	return acos(cosa) * 180 / 3.1415926;
}

double Canvas::calc_angle_sim(int i, int j)
{
	assert(i >= 0 && i < int(poly1.size()));
	assert(j >= 0 && j < int(poly2.size()));

	// triangle1
	double e11 = calc_edge_len(poly1[(poly1.size() + i - 1) % poly1.size()],\
			poly1[i]);
	double e12 = calc_edge_len(poly1[i],\
			poly1[(i + 1) % poly1.size()]);
	double e13 = calc_edge_len(poly1[(poly1.size() + i - 1) % poly1.size()],\
			poly1[(i + 1) % poly1.size()]);

	double a1 = calc_angle(e13, e11, e12);

	// triangle2
	double e21 = calc_edge_len(poly2[(poly2.size() + j - 1) % poly2.size()],\
			poly2[j]);
	double e22 = calc_edge_len(poly2[j],\
			poly2[(j + 1) % poly2.size()]);
	double e23 = calc_edge_len(poly2[(poly2.size() + j - 1) % poly2.size()],\
			poly2[(j + 1) % poly2.size()]);

	double a2 = calc_angle(e23, e21, e22);

	// calculate similarity
	double w1 = 0.5, w2 = 0.5;
	double s = w1 * (1 - abs(e11 * e22 - e12 * e21) / (e11 * e22 + e12 * e21))\
				 + w2 * (1 - abs(a1 - a2) / 360);

	return s;
}

void Canvas::calc_poly_sim()
{
	// size of polygon1 should >= size of polygon2
	assert(poly1.size() >= poly2.size());

	// calculate similarity matrix1
	for(size_t i = 0; i < poly2.size(); ++i) {
		for(size_t j = 0; j < poly1.size(); ++j) {
			sim[i][j] = calc_angle_sim(j, i);
		}
	}
}

void Canvas::calc_points_map()
{
	mapping.clear();
	map<int, int> mapping1;
	double min_dist1 = 10000;
	
	calc_mapping(poly2, mapping1, min_dist1, sim);

	if(min_dist1 < min_dist) {
		min_dist = min_dist1;
		mapping.insert(mapping1.begin(), mapping1.end());
	}
}

void Canvas::calc_mapping(vector<Point> &poly, map<int, int> &_map,\
		double &mdist, double sim_mat[100][100])
{
	for(size_t i = 0; i < poly.size(); ++i) {
		// construct cograph of similarity matrix
		double co_sim_i[100][100];
		for(size_t m = 0; m < poly.size(); ++m) {
			for(size_t n = 0; n < poly1.size(); ++n) {
				co_sim_i[m][n] = 1 - sim_mat[(m + i) % poly.size()][n];
			}
		}
		for(size_t i = 0; i < poly1.size(); ++i) {
			co_sim_i[poly.size()][i] = co_sim_i[0][i];
		}
		for(size_t i = 0; i < poly.size(); ++i) {
			co_sim_i[i][poly1.size()] = co_sim_i[i][0];
		}
		co_sim_i[poly.size()][poly1.size()] = co_sim_i[0][0];

		// calculate shortest path
		vector< pair<int, int> > path;
		double dist = calc_shortest_path(co_sim_i, poly.size() + 1,\
				poly1.size() + 1, path);
		
		// if dist < min_dist, update mapping
		if(dist < mdist) {
			_map.clear();
			for(size_t j = 0; j < path.size(); ++j) {
				pair<int, int> p;
				p.first = path[j].second;
				p.second = (poly.size() + path[j].first - i) % poly.size();
				_map[p.first] = p.second;
			}
			mdist = dist;
		}
	}
}

double Canvas::calc_shortest_path(double mat[100][100], int m, int n,\
		vector< pair<int, int> > &path)
{
	// calculate shortest path
	double dp[100][100];
	dp[0][0] = mat[0][0];

	pair<int, int> prev[100][100];
	
	for(int i = 1; i < m; ++i) {
		dp[i][0] = 100;
		prev[i][0] = make_pair(0, 0);
	}

	for(int i = 1; i < n; ++i) {
		dp[0][i] = dp[0][i - 1] + mat[0][i];
		prev[0][i] = make_pair(0, i - 1);
	}

	for(int i = 1; i < m; ++i) {
		for(int j = 1; j < n; ++j) {
			dp[i][j] = 100;
		}
	}

	for(int i = 1; i < m; ++i) {
		for(int j = 1; j < n; ++j) {
			int x = 0, y = 0;
			if(dp[i - 1][j - 1] <= dp[i][j - 1]) {
				x = i - 1;
				y = j - 1;
			} else {
				x = i;
				y = j - 1;
			}
			dp[i][j] = dp[x][y] + mat[i][j];
			prev[i][j] = make_pair(x, y);
		}
	}

	// get path
	path.clear();
	pair<int, int> p = prev[m - 1][n - 1];
	while(!(p.first == 0 && p.second == 0)) {
		path.push_back(p);
		p = prev[p.first][p.second];
	}
	path.push_back(make_pair(0, 0));
	
	reverse(path.begin(), path.end());

	return dp[m - 1][n - 1];
}

void Canvas::print_mat(double mat[100][100], int m, int n)
{
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			cout << mat[i][j] << " ";
		}
		cout << endl;
	}
}

void Canvas::show_map()
{
	if(show_mapping == true) {
		show_mapping = false;
	} else if(show_mapping == false){
		if(mapping.size() == 0) {
			calc_poly_sim();
			calc_points_map();
			calc_best_affine_trans();
		}
		show_mapping = true;
	}
	update();
}

void Canvas::show_anchor()
{
	if(show_anchor_points == true) {
		show_anchor_points = false;
	} else {
		if(mapping.size() == 0) {
			calc_poly_sim();
			calc_points_map();
		}
		calc_best_affine_trans();
		if(anchor_points.size() != 0) {
			show_anchor_points = true;
		}
	}
	update();
}

double Canvas::calc_angle_smooth(int i)
{
	assert(mapping.size() == poly1.size());
	
	int j = mapping[i];

	vector<Point> angle1, angle2;
	angle1.push_back(poly1[(poly1.size() + i - 1) % poly1.size()]);
	angle1.push_back(poly1[i]);
	angle1.push_back(poly1[(i + 1) % poly1.size()]);
	
	angle2.push_back(poly2[(poly2.size() + j - 1) % poly2.size()]);
	angle2.push_back(poly2[j]);
	angle2.push_back(poly2[(j + 1) % poly2.size()]);

	// triangle1
	double e11 = calc_edge_len(angle1[0], angle1[1]);
	double e12 = calc_edge_len(angle1[1], angle1[2]);
	double e13 = calc_edge_len(angle1[0], angle1[2]);

	double a11 = calc_angle(e13, e11, e12);
	double a12 = calc_angle(e11, e13, e12);
	double a13 = calc_angle(e12, e13, e11);

	// triangle2
	double e21 = calc_edge_len(angle2[0], angle2[1]);
	double e22 = calc_edge_len(angle2[1], angle2[2]);
	double e23 = calc_edge_len(angle2[0], angle2[2]);

	double a21 = calc_angle(e23, e21, e22);
	double a22 = calc_angle(e21, e23, e22);
	double a23 = calc_angle(e22, e23, e21);
	

	// smooth = a * S + b * (1 - R / 180) + c * A
	double smooth = 0;
	double a = 0.4, b = 0.3, c = 0.3;

	// calculate S
	double w1 = 0.5, w2 = 0.5;
	double S = w1 * (1 - (abs(e11 - e21) + abs(e12 - e22) + \
				abs(e13 - e23)) / (e11 + e21 + e12 + e22 + e13 + e23))\
				+ w2 * (1 - (abs(a11 - a21) + abs(a12 - a22) + abs(a13 - a23)) / 180);
	
	// calculate R
	Eigen::Matrix2d A;
	Eigen::Vector2d T;
	calc_affine_mat(angle1, angle2, A, T);
	Eigen::Matrix2d B, C;
	decompose_affine_mat(A, B, C);
	double R = acos(B(0, 0)) * 180 / 3.14159262;

	// calculate A
	double angle_area1 = calc_area(angle1);
	double angle_area2 = calc_area(angle2);
	double poly_area1 = calc_area(poly1);
	double poly_area2 = calc_area(poly2);
	double _a = (angle_area1 + angle_area2) / (poly_area1 + poly_area2);

	smooth = a * S + b * (1 - R / 180) + c * _a;

	if(a11 > 180 || a21 > 180) {
		smooth = 0;
	}
	return smooth;
}

double Canvas::calc_area(vector<Point> &poly)
{
	double area = 0;
	for(size_t i = 0; i < poly.size(); ++i) {
		area += 0.5 * (poly[i].x * poly[(i + 1) % poly.size()].y\
				- poly[(i + 1) % poly.size()].x * poly[i].y);
	}
	return fabs(area);
}

double Canvas::calc_smooth(int a, int b, int c)
{
	return calc_angle_smooth(a) * calc_angle_smooth(b) * calc_angle_smooth(c);	
}

void Canvas::calc_best_affine_trans()
{
	assert(mapping.size() == poly1.size());

	vector<int> perm;
	for(size_t i = 0; i < mapping.size() - 3; ++i){
		perm.push_back(0);
	}
	perm.push_back(1);
	perm.push_back(1);
	perm.push_back(1);

	double max_smooth = 0;

	bool flag = true;
	while(flag) {
		vector<int> points;
		vector<int>::iterator it = perm.begin() - 1;
		for(int i = 0; i < 3; ++i) {
			it = ::find(it + 1, perm.end(), 1);
			points.push_back(it - perm.begin());
		}

		flag = next_permutation(perm.begin(), perm.end());

		// check if 2 points mapping to the same point in polygon2
		if(mapping[points[0]] == mapping[points[1]] ||\
				mapping[points[0]] == mapping[points[2]] ||\
				mapping[points[1]] == mapping[points[2]]) {
			continue;
		}

		double smooth = calc_smooth(points[0], points[1], points[2]);
		if(smooth > max_smooth) {
			anchor_points.clear();
			anchor_points.assign(points.begin(), points.end());
			max_smooth = smooth;
		}
	}
	
	// calculate affine matrix
	vector<Point> angle1, angle2;
	for(auto i : anchor_points) {
		angle1.push_back(poly1[i]);
		angle2.push_back(poly2[mapping[i]]);
	}

	calc_affine_mat(angle1, angle2, _A, _T);
	decompose_affine_mat(_A, _B, _C);
}

void Canvas::calc_affine_mat(vector<Point> &angle1, vector<Point> &angle2,\
		Eigen::Matrix2d &A, Eigen::Vector2d &T)
{
	Eigen::Matrix<double, 3, 3> before, after;

	before << angle1[0].x, angle1[1].x, angle1[2].x,\
			  angle1[0].y, angle1[1].y, angle1[2].y,\
			  1,          1,          1;

	after << angle2[0].x, angle2[1].x, angle2[2].x,\
			 angle2[0].y, angle2[1].y, angle2[2].y,\
			 1,          1,          1;
	
	// transform parameters matrix
	Eigen::Matrix<double, 3, 3> paras;
	paras = after * before.inverse();

	// composite matrix
	A << paras(0, 0), paras(1, 0),\
		 paras(0, 1), paras(1, 1);
	// translation matrix
	T << paras(0, 2), paras(1, 2);
}

void Canvas::decompose_affine_mat(Eigen::Matrix2d &A, Eigen::Matrix2d &B,\
		Eigen::Matrix2d &C)
{
	// B = A + sign(det(A)) * [a22, -a21
	// 						   -a12, a11]
	// C = B^-1 * A
	
	double det_A = A.determinant();
	double sign_det_A = det_A > 0 ? 1 : -1;
	if(det_A == 0) {
		sign_det_A = 0;
	}

	Eigen::Matrix2d temp;
	temp << A(1, 1), -A(1, 0),\
			-A(0, 1), A(0, 0);

	B = A + sign_det_A * temp;
	C = B.inverse() * A;
	double t = sqrt(B(0, 0) * B(0, 0) + B(0, 1) * B(0, 1));
	B /= t;
	C *= t;
}

vector< vector<Point> > Canvas::interpolation()
{
	assert(mapping.size() == poly1.size());

	vector< vector<Point> > frames;
	vector<Point> origin_local, dest_local;
	vector<Point> origin_anchor, dest_anchor;
	for(auto i : anchor_points) {
		origin_anchor.push_back(poly1[i]);
		dest_anchor.push_back(poly2[mapping[i]]);
	}

	for(auto i : mapping) {
		origin_local.push_back(calc_local_coords(origin_anchor, poly1[i.first]));
		dest_local.push_back(calc_local_coords(dest_anchor, poly2[i.second]));
	}

	for(double t = 0; t <= 1; t += 0.01) {
		vector<Point> frame;
		for(size_t i = 0; i < origin_local.size(); ++i) {
			Point p;
			double x = (1 - t) * origin_local[i].x + t * dest_local[i].x;
			double y = (1 - t) * origin_local[i].y + t * dest_local[i].y;
			
			vector<Point> anchor = interpolate_anchor(origin_anchor, t);
			p.x = anchor[1].x + x * (anchor[0].x - anchor[1].x) + \
				  y * (anchor[2].x - anchor[1].x);
			p.y = anchor[1].y + x * (anchor[0].y - anchor[1].y) + \
				  y * (anchor[2].y - anchor[1].y);
			frame.push_back(p);
		}
		frames.push_back(frame);
	}

	return frames;
}

vector<Point> Canvas::interpolate_anchor(vector<Point> &origin_anchor, double t)
{
	Eigen::Matrix2d temp;
	temp << 1, 0,\
			0, 1;

	Eigen::Matrix2d rotate;
	rotate << cos(t * acos(_B(0, 0))), -sin(t * asin(-_B(0, 1))),\
		 sin(t * asin(_B(1, 0))), cos(t * acos(_B(1, 1)));

	Eigen::Matrix2d scale;
	scale << t * _C(0, 0), t * _C(0, 1),\
			 t * _C(1, 0), t * _C(1, 1);
	
	Eigen::Vector2d trans;
	trans << t * _T(0), t * _T(1);

	Eigen::Matrix2d a = (1 - t) * temp + rotate * scale;

	vector<Point> anchor;
	for(size_t i = 0; i < origin_anchor.size(); ++i) {
		Eigen::Vector2d coord1;
		coord1 << origin_anchor[i].x, origin_anchor[i].y;
		Eigen::Vector2d p = coord1.transpose() * a + trans.transpose();
		anchor.push_back(Point(p(0), p(1)));
	}

	return anchor;
}

Point Canvas::calc_local_coords(vector<Point> &angle, Point &p)
{
	Point coord;
	
	// solve equation set
	Eigen::Matrix3d trans;
	trans << angle[0].x - angle[1].x, angle[2].x - angle[1].x, angle[1].x,\
			 angle[0].y - angle[1].y, angle[2].y - angle[1].y, angle[1].y,\
			 0, 0, 1;

	Eigen::Vector3d after;
	after << p.x, p.y, 1;

	Eigen::Vector3d before = trans.inverse() * after;
	coord.x = before(0);
	coord.y = before(1);

	return coord;
}
