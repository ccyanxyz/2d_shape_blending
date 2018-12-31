#include <QApplication>
#include <QObject>
#include <QWidget>
#include <QSpinBox>
#include <QSlider>
#include <QHBoxLayout>
#include <QGridLayout>
#include <QLabel>
#include <QPushButton>
#include <QPainter>
#include <QLineEdit>
#include <QDebug>
#include <QMessageBox>
#include "canvas.h"

#include <iostream>

class Window: public QWidget
{
	Q_OBJECT
	private:
		Canvas *canvas; // canvas to draw spline and car
		QWidget *window; // menu window
		QVBoxLayout *layout;
		QHBoxLayout *hlayout;

		QPushButton *clear_button;
		QPushButton *play_button;
		QPushButton *draw_done;

	public:
		Window(QWidget *parent = 0): QWidget(parent)
		{
			canvas = new Canvas();
			QWidget *win = new QWidget();

			draw_done = new QPushButton(this);
			draw_done->setText("Done");
			clear_button = new QPushButton(this);
			clear_button->setText("Clear");
			play_button = new QPushButton(this);
			play_button->setText("Play");

			connect(draw_done, SIGNAL(clicked()), this, SLOT(draw_end()));
			connect(clear_button, SIGNAL(clicked()), this, SLOT(clear_canvas()));
			connect(play_button, SIGNAL(clicked()), this, SLOT(play()));

			layout = new QVBoxLayout();
			hlayout = new QHBoxLayout();
			layout->addWidget(draw_done);
			layout->addWidget(clear_button);
			layout->addWidget(play_button);
			
			win->setFixedWidth(204);
			win->setLayout(layout);
			this->resize(800, 600);
			hlayout->addWidget(win);
			hlayout->addWidget(canvas);
			this->setLayout(hlayout);
		}

	private slots:
		// finish 1 draw
		void draw_end()
		{
			canvas->draw_done();
		}

		// clear canvas when clear button is clicked
		void clear_canvas()
		{
			canvas->clear();
			update();
		}

		// start to move the car when play button is clicked
		void play()
		{
			canvas->play();
		}
};
