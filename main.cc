#include <QApplication>
#include "window.hpp"

int main(int argc, char *argv[])
{
	QApplication app(argc, argv);

	Window win;
	win.show();

	app.exec();

	return 0;
}
