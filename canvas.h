#include <QWidget>
#include <map>
#include <Eigen/Dense>

using namespace std;
class Point {

public:
	int x;
	int y;

	Point(int _x, int _y): x(_x), y(_y) {  }
	Point() {
		x = 0;
		y = 0;
	}
};

class Canvas: public QWidget
{
	Q_OBJECT
	private:
		vector<Point> poly1;
		vector<Point> poly2;
		int num;
		// similarity matrix
		double sim[100][100];
		// shortest path length
		double min_dist;
		// shortest path points mapping
		//vector< pair<int, int> > mapping;
		map<int, int> mapping;
		bool show_mapping;

		// anchor points
		vector<int> anchor_points;
		// show anchor
		bool show_anchor_points;

	public:
		// contrustor
		Canvas(QWidget *parent = 0): QWidget(parent)
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
		}

		void clear();

		void draw_done();

		void play();

		void show_map();

		// calculate edge length
		double calc_edge_len(Point &a, Point &b);
		// calculate angle value
		double calc_angle(double a, double b, double c);
		// calculate angle similarity
		double calc_angle_sim(int i, int j);
		// calculate polygon similarity
		void calc_poly_sim();
		// calculate mapping
		void calc_points_map();
		// calculate shortest path
		double calc_shortest_path(double mat[100][100], int m, int n,\
				vector< pair<int, int> > &path);
		// print matrix
		void print_mat(double mat[100][100], int m, int n);

		// calculate smooth of angle pair (i, mapping[i])
		double calc_angle_smooth(int i);

		// calculate smooth of 3 paired angles
		double calc_smooth(int a, int b, int c);
		// calculate polygon area
		double calc_area(vector<Point> &poly);

		// get 3 pair of angles with biggest smooth value
		// then calculate matrix A, B, C
		void calc_best_affine_trans(Eigen::Matrix2d &, Eigen::Matrix2d &,
				Eigen::Matrix2d &, Eigen::Vector2d &);

		// calculate affine transform matrix
		void calc_affine_mat(vector<Point> &, vector<Point> &, Eigen::Matrix2d &,\
				Eigen::Vector2d &);

		// decompose affine transform matrix
		void decompose_affine_mat(Eigen::Matrix2d &, Eigen::Matrix2d &,\
				Eigen::Matrix2d &);

		// show anchor points
		void show_anchor();
		// interpolation
		void interpolation();

	protected:
		void paintEvent(QPaintEvent *);
			
		void mousePressEvent(QMouseEvent *e);
};
