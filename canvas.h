#include <QWidget>

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
		vector< pair<int, int> > mapping;
		bool show_mapping;

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

	protected:
		void paintEvent(QPaintEvent *);
			
		void mousePressEvent(QMouseEvent *e);
};
