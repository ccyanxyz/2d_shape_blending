#include <iostream>
#include <algorithm>
#include <QWidget>
#include <QTimer>
#include <QPainter>
#include <QBrush>
#include <QPixmap>
#include <QMouseEvent>
#include <QDebug>
#include <QMessageBox>
#include <QLabel>
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
}

void Canvas::draw_done()
{
	if(num <= 1) {
		num += 1;
		update();
	}
}

void Canvas::play()
{
	calc_poly_sim();
	calc_points_map();
}

void Canvas::paintEvent(QPaintEvent *)
{
	QPainter p(this);

	QBrush brush(Qt::black, Qt::SolidPattern);
	p.setPen(QColor(Qt::red));
	p.setBrush(brush);
	// paint control points
	for(size_t i = 0; i < poly1.size(); ++i){
		p.drawEllipse(poly1[i].x, poly1[i].y, 5, 5);
	}

	// paint line
	for(size_t i = 0; i < poly1.size(); ++i){
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
	if(show_mapping && mapping.size() == poly1.size()) {
		p.setPen(QColor(Qt::green));
		for(size_t i = 0; i < mapping.size(); ++i) {
			p.drawLine(poly1[mapping[i].first].x, poly1[mapping[i].first].y,\
					poly2[mapping[i].second].x, poly2[mapping[i].second].y);
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

	// calculate similarity matrix
	for(size_t i = 0; i < poly2.size(); ++i) {
		for(size_t j = 0; j < poly1.size(); ++j) {
			sim[i][j] = calc_angle_sim(j, i);
		}
	}

	cout << "sim:" << endl;
	print_mat(sim, poly2.size(), poly1.size());
}

void Canvas::calc_points_map()
{
	mapping.clear();
	for(size_t i = 0; i < poly2.size(); ++i) {
		// construct cograph of similarity matrix
		double co_sim_i[100][100];
		for(size_t m = 0; m < poly2.size(); ++m) {
			for(size_t n = 0; n < poly1.size(); ++n) {
				co_sim_i[m][n] = 1 - sim[(m + i) % poly2.size()][n];
			}
		}
		for(size_t i = 0; i < poly1.size(); ++i) {
			co_sim_i[poly2.size()][i] = co_sim_i[0][i];
		}
		for(size_t i = 0; i < poly2.size(); ++i) {
			co_sim_i[i][poly1.size()] = co_sim_i[i][0];
		}
		co_sim_i[poly2.size()][poly1.size()] = co_sim_i[0][0];

		//cout << "co_sim:" << endl;
		//print_mat(co_sim_i, poly2.size() + 1, poly1.size() + 1);

		// calculate shortest path
		vector< pair<int, int> > path;
		double dist = calc_shortest_path(co_sim_i, poly2.size() + 1,\
				poly1.size() + 1, path);
		
		// if dist < min_dist, update mapping
		if(dist < min_dist) {
			mapping.clear();
			for(size_t j = 0; j < path.size(); ++j) {
				pair<int, int> p;
				p.first = path[j].second;
				p.second = (poly2.size() + path[j].first - i) % poly2.size();
				mapping.push_back(p);
			}
			min_dist = dist;
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

	/*
	cout << "dp:" << endl;
	print_mat(dp, m, n);

	cout << "prev:" << endl;
	for(int i = 0; i < m; ++i) {
		for(int j = 0; j < n; ++j) {
			cout << "(" << prev[i][j].first << ", " << prev[i][j].second << ") ";
		}
		cout << endl;
	}
	cout << endl;
	*/

	// get path
	path.clear();
	pair<int, int> p = prev[m - 1][n - 1];
	while(p.first != 0 && p.second != 0) {
		path.push_back(p);
		p = prev[p.first][p.second];
	}
	path.push_back(make_pair(0, 0));
	
	reverse(path.begin(), path.end());

	for(size_t i = 0; i < path.size(); ++i) {
		cout << "(" << path[i].first << ", " << path[i].second << ") ";
	}
	cout << dp[m-1][n-1];
	cout << endl;

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
		calc_poly_sim();
		calc_points_map();
		show_mapping = true;
	}
	cout << "===== mapping =====" << endl;
	for(size_t i = 0; i < mapping.size(); ++i) {
		cout << mapping[i].first << "->" << mapping[i].second << endl;
	}
	cout << "===================" << endl;
	update();
}
