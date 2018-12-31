#include <iostream>
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
	min_dist = 100;
	path.clear();
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
	print_mat(sim, poly2.size() + 1, poly1.size() + 1);
}

void Canvas::calc_points_map()
{
	// cograph of similarity matrix
	double co_sim[100][100];
	for(size_t i = 0; i < poly2.size() + 1; ++i) {
		for(size_t j = 0; j < poly1.size() + 1; ++j) {
			co_sim[i][j] = 1 - sim[i][j];
		}
	}

	cout << "sim:" << endl;
	print_mat(sim, poly2.size() + 1, poly1.size() + 1);

	cout << "co_graph:" << endl;
	print_mat(co_sim, poly2.size() + 1, poly1.size() + 1);
	// calculate shortest path of co_sim
	calc_shortest_path(co_sim, poly2.size() + 1, poly1.size() + 1);	
}

double Canvas::calc_shortest_path(double mat[100][100], int m, int n)
{
	double dp[100][100];
	dp[0][0] = mat[0][0];
	for(int i = 1; i < m; ++i) {
		dp[i][0] = 100;
		dp[0][i] = dp[0][i - 1] + mat[0][i];
	}

	for(int i = 1; i < m; ++i) {
		for(int j = 1; j < m; ++j){
			dp[i][j] = 100;
		}
	}

	vector< pair<int, int> > path1;

	for(int i = 1; i < m; ++i) {
		for(int j = 1; j < n; ++j) {
			int x = 0, y = 0;
			if(mat[i - 1][j - 1] <= mat[i - 1][j]) {
				x = i - 1;
				y = j - 1;
			} else if(mat[i - 1][j] <= mat[i - 1][j - 1]) {
				x = i - 1;
				y = j;
			}

			if(mat[x][y] + mat[i][j] < dp[i][j]) {
				dp[i][j] = mat[x][y] + mat[i][j];
				path1.push_back(make_pair(x, y));
			}
		}
	}

	cout << "dp:" << endl;
	print_mat(dp, m, n);

	cout << "path:" << endl;
	for(size_t i = 0; i < path1.size(); ++i) {
		cout << path1[i].first << ", " << path1[i].second << endl;
	}

	return dp[m][n];
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
