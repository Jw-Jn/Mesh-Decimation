/*
	Based on examples of NanoGUI.
	NanoGUI was developed by Wenzel Jakob <wenzel.jakob@epfl.ch>.
	The widget drawing code is based on the NanoVG demo application
	by Mikko Mononen.
	All rights reserved. Use of this source code is governed by a
	BSD-style license that can be found in the LICENSE.txt file.
*/

////// add while windows//////////
// #define _CRT_SECURE_NO_DEPRECATE

// #include <GL/glew.h>
///////////////////////////////

#include <nanogui/opengl.h>
#include <nanogui/glutil.h>
#include <nanogui/screen.h>
#include <nanogui/window.h>
#include <nanogui/layout.h>
#include <nanogui/label.h>
#include <nanogui/checkbox.h>
#include <nanogui/button.h>
#include <nanogui/toolbutton.h>
#include <nanogui/popupbutton.h>
#include <nanogui/combobox.h>
#include <nanogui/progressbar.h>
#include <nanogui/entypo.h>
#include <nanogui/messagedialog.h>
#include <nanogui/textbox.h>
#include <nanogui/slider.h>
#include <nanogui/imagepanel.h>
#include <nanogui/imageview.h>
#include <nanogui/vscrollpanel.h>
#include <nanogui/colorwheel.h>
#include <nanogui/graph.h>
#include <nanogui/tabwidget.h>
#include <nanogui/glcanvas.h>
#include <iostream>
#include <string>

#include <cstdint>
#include <memory>
#include <utility>

#include "WingEdge.h"


#if defined(__GNUC__)
#  pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#endif
#if defined(_WIN32)
#  pragma warning(push)
#  pragma warning(disable: 4457 4456 4005 4312)
#endif

#define STB_IMAGE_IMPLEMENTATION
#include <stb_image.h>

#if defined(_WIN32)
#  pragma warning(pop)
#endif
#if defined(_WIN32)
#  if defined(APIENTRY)
#    undef APIENTRY
#  endif
#  include <windows.h>
#endif

using std::cout;
using std::cerr;
using std::endl;
using std::string;
using std::vector;
using std::pair;
using std::to_string;

using nanogui::Screen;
using nanogui::Window;
using nanogui::GroupLayout;
using nanogui::Button;
using nanogui::CheckBox;
using nanogui::Vector2f;
using nanogui::Vector2i;
using nanogui::MatrixXu;
using nanogui::MatrixXf;
using nanogui::Label;
using nanogui::Arcball;

class MyGLCanvas : public nanogui::GLCanvas {
public:
	WingEdge mesh;

	MyGLCanvas(Widget *parent) : nanogui::GLCanvas(parent) {
		using namespace nanogui;

		mShader.initFromFiles("a_smooth_shader", "StandardShading.vertexshader", "StandardShading.fragmentshader");
		
		mesh.loadOBJfile("./obj/cube.obj");
		updateMesh();

		mShader.bind();

		mShader.uploadAttrib("vertexPosition_modelspace", positions);
		mShader.uploadAttrib("color", colors);
		mShader.uploadAttrib("vertexNormal_modelspace", normals);

		Matrix4f V;
		V.setIdentity();
		mShader.setUniform("V", V);

		Matrix4f M;
		M.setIdentity();
		mShader.setUniform("M", M);

		mShader.setUniform("LightPosition_worldspace", Vector3f(-2, 6, -4));

		mScale.setIdentity();
		mRotation.setIdentity();
		mTranslation.setIdentity();
		renderMode = 0;
		collapseNum = 1;
		k = 1;
		transX = 0;
		transY = 0;
		rotX = 0;
		rotY = 0;
		rotZ = 0;
		//divMethod = "Butterfly";
	}

	~MyGLCanvas() {
		mShader.free();
	}

	void updateMeshPositions(MatrixXf newPositions) {
		positions = newPositions;
	}

	void setRotation() {
		mRotation.setIdentity();
		mRotation.topLeftCorner<3, 3>() = Eigen::Matrix3f(Eigen::AngleAxisf(rotX, nanogui::Vector3f::UnitX()) *
			Eigen::AngleAxisf(rotY, nanogui::Vector3f::UnitY()) *
			Eigen::AngleAxisf(rotZ, nanogui::Vector3f::UnitZ()));
	}

	void setRotationX(float x) {
		rotX = x;
	}

	void setRotationY(float y) {
		rotY = y;
	}

	void setRotationZ(float z) {
		rotZ = z;
	}

	void setRenderMode(int a){
		renderMode = a;
	}

	void setScale(float vScale) {
		mScale.setIdentity();
		mScale.topLeftCorner<3, 3>() *= vScale;
	}

	void setTranslate() {
		mTranslation.setIdentity();
		mTranslation.topRightCorner<3, 1>() = nanogui::Vector3f(transX, transY, 0.0f);
	}

	void setTranslateX(float x) {
		transX = x;
	}

	void setTranslateY(float y) {
		transY = y;
	}

	/*void setDivMethod(int x) {
		std::cout << "method" << x << std::endl;
		if (x == 0)
			divMethod = "Butterfly";
		else
			divMethod = "Loop";
	}*/

	/*void subdivision() {
		mesh.subdivision(divMethod);
	}*/

	void setMultipleChoice(int a){
		k = a;
	}

	void setCollapseNum(int a){
		collapseNum = a;
	}

	void decimation() {

		if ((collapseNum*2) < faceNum)
		{
			bool result = mesh.decimation(k, collapseNum);
			if (!result)
				std::cout<<"nothing to collapse!"<<std::endl;
		}
		else
			std::cout<<"nothing to collapse!"<<std::endl;
	}

	void updateMesh(){
		std::vector<std::vector<Wvertex *>> faces = mesh.extractVerticesOfFaces();
		faceNum = faces.size();

		std::cout << "faceNum: " <<faceNum << std::endl;

		positions = MatrixXf(3, faceNum *9);
		normals = MatrixXf(3, faceNum *9);
		colors = MatrixXf(3, faceNum *9);

		int k = 0;
		for (int i = 0; i < faceNum; ++i) {
			for (int j = 0; j < faces[i].size(); ++j) {
				vec3 v = mesh.getVertex(faces[i][j]);
				vec3 n = mesh.getNorm(faces[i][j]);
				positions.col(k) << v.x, v.y, v.z;
				normals.col(k) << n.x, n.y, n.z;
				colors.col(k) << 1.0, 0, 0;
				k++;
			}
		}

		for (int i = 0; i < faceNum; ++i) {
			for (int j = 0; j < faces[i].size(); ++j) {
				vec3 n = mesh.getNorm(faces[i][j]);
				vec3 v = mesh.getVertex(faces[i][j]) + n * 0.001;
				positions.col(k) << v.x, v.y, v.z;
				normals.col(k) << n.x, n.y, n.z;
				colors.col(k) << 0, 0, 0;
				k++;

				n = mesh.getNorm(faces[i][(j + 1) % faces[i].size()]);
				v = mesh.getVertex(faces[i][(j+1)% faces[i].size()]) + n * 0.001;
				positions.col(k) << v.x, v.y, v.z;
				normals.col(k) << n.x, n.y, n.z;
				colors.col(k) << 0, 0, 0;
				k++;
			}
		}
	}

	virtual void drawGL() override {
		using namespace nanogui;

		mShader.bind();

		mShader.uploadAttrib("vertexPosition_modelspace", positions);
		mShader.uploadAttrib("color", colors);
		mShader.uploadAttrib("vertexNormal_modelspace", normals);

		Matrix4f mvp;
		mvp.setIdentity();
		mvp.topLeftCorner<3, 3>() = mvp.topLeftCorner<3, 3>();
		setTranslate();
		setRotation();
		mvp = mTranslation * mRotation * mScale * mvp;
		mShader.setUniform("MVP", mvp);

		Matrix4f V;
		V.setIdentity();
		mShader.setUniform("V", V);

		Matrix4f M;
		M.setIdentity();
		mShader.setUniform("M", M);

		glEnable(GL_DEPTH_TEST);

		MatrixXf flat_normals;

		switch (renderMode)
		{	
			// with edge
			case 0:
				mShader.drawArray(GL_TRIANGLES, 0, faceNum * 3);
				mShader.drawArray(GL_LINES, faceNum * 3, faceNum * 9);
			break;
			// wireframe
			case 1:
				mShader.drawArray(GL_LINES, faceNum * 3, faceNum * 9);
			break;
			
			// flat shaded
			case 2:	
				flat_normals = MatrixXf(3, faceNum * 3);
				for (int i = 0; i < faceNum * 3; i++)
				{
					if (i % 3 == 0)
						flat_normals.col(i) = (normals.col(i) + normals.col(i + 1) + normals.col(i + 2)) / 3;
					else
						flat_normals.col(i) = flat_normals.col(i - 1);
				}
				mShader.uploadAttrib("vertexNormal_modelspace", flat_normals);

				mShader.drawArray(GL_TRIANGLES, 0, faceNum * 3);
				break;

			// smooth shaded
			case 3:	
				mShader.drawArray(GL_TRIANGLES, 0, faceNum * 3);
				break;

		}

		glDisable(GL_DEPTH_TEST);
	}

	
private:
	MatrixXf positions;
	MatrixXf normals;
	MatrixXf colors;

	nanogui::GLShader mShader;
	Eigen::Matrix4f mRotation;
	Eigen::Matrix4f mScale;
	Eigen::Matrix4f mTranslation;
	float transX, transY, rotX, rotY, rotZ;
	int renderMode;
	int faceNum;
	int k;
	int collapseNum;
	//std::string divMethod;
};


class ExampleApplication : public nanogui::Screen {
public:
	ExampleApplication() : nanogui::Screen(Eigen::Vector2i(800, 530), "Mesh Subdivision", false) {
		using namespace nanogui;
		
		Window *window = new Window(this, "");
		window->setPosition(Vector2i(320, 40));
		// window->setFixedSize(Vector2i(450, 450));
		window->setLayout(new GroupLayout());

		///////////////// add while windows ///////////////////
		// glewExperimental = GL_TRUE;
		// GLenum err = glewInit();

		// if (GLEW_OK != err)
		// {
		// 	/* Problem: glewInit failed, something is seriously wrong. */
		// 	fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
		// }
		///////////////////////////////////////////////////


		mCanvas = new MyGLCanvas(window);
		mCanvas->setBackgroundColor({ 100, 100, 100, 255 });
		mCanvas->setSize({ 400, 400 });

		Widget *tools = new Widget(window);
		tools->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 5));

		nanogui::GLShader mShader;

		Window *anotherWindow = new Window(this, "");
		anotherWindow->setPosition(Vector2i(30, 15));
		anotherWindow->setLayout(new GroupLayout());

		new Label(anotherWindow, "File", "sans-bold");
		tools = new Widget(anotherWindow);
		tools->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 6));
		Button *b = new Button(tools, "Open");
		b->setCallback([&] {
			bool result = mCanvas->mesh.loadOBJfile(file_dialog({ {"obj", "obj file"} }, false));
			if (result){
				cout << "File loaded!" << endl;
				mCanvas->updateMesh();
			}
			//cout << "File dialog result: " << file_dialog({ {"obj", "obj file"} }, false) << endl;
		});

		b = new Button(tools, "Save");
		b->setCallback([&] {
			bool result = mCanvas->mesh.saveToOBJfile(file_dialog({ {"obj", "obj file"} }, true));
			if (result){
				cout << "File saved!" << endl;
				auto dlg = new MessageDialog(this, MessageDialog::Type::Information, "Info", "File saved!");
				dlg->setCallback([](int result) { cout << "Dialog result: " << result << endl; });
			}
			// cout << "File dialog result: " << file_dialog({ {"obj", "obj file"} }, true) << endl;
		});

		new Label(anotherWindow, "Render options", "sans-bold");
		ComboBox *combo = new ComboBox(anotherWindow, {"shaded with mesh edges", "wireframe", "flat shaded", "smooth shaded"});
		combo->setCallback([&](int value) {
			mCanvas->setRenderMode(value);
		});

		new Label(anotherWindow, "Rotate", "sans-bold");

		Widget *panelRot = new Widget(anotherWindow);
		panelRot->setLayout(new GridLayout(Orientation::Horizontal, 2,
			Alignment::Middle, 0, 0));

		{new Label(panelRot, "x :", "sans-bold");
		Slider *rotSliderX = new Slider(panelRot);
		rotSliderX->setValue(0.5f);
		rotSliderX->setFixedWidth(150);
		rotSliderX->setCallback([&](float value) {
			float radians = (value - 0.5f) * 2 * 2 * 3.1415;
			mCanvas->setRotationX(radians);
		}); }

		{new Label(panelRot, "y :", "sans-bold");
		Slider *rotSliderY = new Slider(panelRot);
		rotSliderY->setValue(0.5f);
		rotSliderY->setFixedWidth(150);
		rotSliderY->setCallback([&](float value) {
			float radians = (value - 0.5f) * 2 * 2 * 3.1415;
			mCanvas->setRotationY(radians);
		}); }

		{new Label(panelRot, "z :", "sans-bold");
		Slider *rotSliderZ = new Slider(panelRot);
		rotSliderZ->setValue(0.5f);
		rotSliderZ->setFixedWidth(150);
		rotSliderZ->setCallback([&](float value) {
			float radians = (value - 0.5f) * 2 * 2 * 3.1415;
			mCanvas->setRotationZ(radians);
		}); }

		// translation slider
		new Label(anotherWindow, "Translate", "sans-bold");

		Widget *panelTrans = new Widget(anotherWindow);
		panelTrans->setLayout(new GridLayout(Orientation::Horizontal, 2,
			Alignment::Middle, 0, 0));

		{new Label(panelTrans, "x :", "sans-bold");
		Slider *transSliderX = new Slider(panelTrans);
		transSliderX->setValue(0.5f);
		transSliderX->setFixedWidth(150);
		transSliderX->setCallback([&](float value) {
			mCanvas->setTranslateX((value - 0.5f) * 2);
		}); }

		{new Label(panelTrans, "y :", "sans-bold");
		Slider *transSliderY = new Slider(panelTrans);
		transSliderY->setValue(0.5f);
		transSliderY->setFixedWidth(150);
		transSliderY->setCallback([&](float value) {
			mCanvas->setTranslateY((value - 0.5f) * 2);
		}); }

		new Label(anotherWindow, "Zoom", "sans-bold");

		Widget *panel = new Widget(anotherWindow);
		panel->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 20));

		Slider *slider = new Slider(panel);
		slider->setValue(0.5f);
		slider->setFixedWidth(100);
		TextBox *textBox = new TextBox(panel);
		textBox->setFixedSize(Vector2i(60, 25));
		textBox->setValue("100");
		textBox->setUnits("%");
		slider->setCallback([textBox](float value) {
			textBox->setValue(std::to_string((int)((value*2) * 100)));
		});
		slider->setFinalCallback([&](float value) {
			mCanvas->setScale(value*2);
		});
		textBox->setFontSize(20);
		textBox->setAlignment(TextBox::Alignment::Right);

		// mesh decimation
		new Label(anotherWindow, "Decimation", "sans-bold");

		Widget *panelDec = new Widget(anotherWindow);
		panelDec->setLayout(new BoxLayout(Orientation::Horizontal,
			Alignment::Middle, 0, 10));

		{new Label(panelDec, "k :", "sans-bold");
		TextBox *textBoxK = new TextBox(panelDec);
		textBoxK->setFixedSize(Vector2i(45, 25));
		textBoxK->setEditable(true);
		textBoxK->setValue("1");
		textBoxK->setFontSize(20);
		textBoxK->setAlignment(TextBox::Alignment::Right);
		textBoxK->setCallback([&](const std::string &value) {
			std::stringstream ss(value); 
			int i;
			ss >> i;
			mCanvas->setMultipleChoice(i);
            // std::cout << "Current Manual Value: " << value << std::endl;
            // if (value == "no")
            //     return false;
            return true;
        });}

		{new Label(panelDec, "collapse num :", "sans-bold");
		TextBox *textBoxEdges = new TextBox(panelDec);
		textBoxEdges->setFixedSize(Vector2i(45, 25));
		textBoxEdges->setEditable(true);
		textBoxEdges->setValue("1");
		textBoxEdges->setFontSize(20);
		textBoxEdges->setAlignment(TextBox::Alignment::Right);
		textBoxEdges->setCallback([&](const std::string &value) {
			std::stringstream ss(value); 
			int i;
			ss >> i;
			mCanvas->setCollapseNum(i);
            // std::cout << "Current c Value: " << value << std::endl;
            // if (value == "no")
            //     return false;
            return true;
        });}

		b = new Button(anotherWindow, "Start Decimation");
		b->setCallback([&] {
			mCanvas->decimation();
			mCanvas->updateMesh();
		});

		// QUIT
		b = new Button(anotherWindow, "Quit");
		b->setCallback([&] {
			shutdown();
		});
		
		performLayout();
	}

	virtual void drawContents() override {
	}

	virtual void draw(NVGcontext *ctx) {
		Screen::draw(ctx);
	}


private:
	nanogui::ProgressBar *mProgress;
	MyGLCanvas *mCanvas;
	nanogui::Vector3f rotate_start, rotate_end;
	nanogui::Vector3f translate_start, translate_end;
	int flag_r, flag_t;
};

int main(int /* argc */, char ** /* argv */) {
	try {
		nanogui::init();

		/* scoped variables */ {
			nanogui::ref<ExampleApplication> app = new ExampleApplication();
			app->drawAll();
			app->setVisible(true);
			nanogui::mainloop();
		}

		nanogui::shutdown();
	}
	catch (const std::runtime_error &e) {
		std::string error_msg = std::string("Caught a fatal error: ") + std::string(e.what());
	#if defined(_WIN32)
		MessageBoxA(nullptr, error_msg.c_str(), NULL, MB_ICONERROR | MB_OK);
	#else
		std::cerr << error_msg << endl;
	#endif
		return -1;
	}

	return 0;
}
