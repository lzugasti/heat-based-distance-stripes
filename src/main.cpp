// LEOPOLDO ZUGASTI 260919951

#pragma once
#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <string>

#include "GLHeaders.h"

#include "GLSL.h"
#include "MatrixStack.h"
#include "Program.h"

#include "Face.h"
#include "HalfEdge.h"
#include "HEDS.h"
#include "PolygonSoup.h"
#include "MeshDrawHeatGeo.h"
#include "MeshDrawPicking.h"
#include "ArcBall.h"

using namespace std;

/* heat solution time, feel free to adjust */
double tSol = 1e2; // [1e-3, 1e3]
/* heat solution step, feel free to adjust */
int steps = 1; // [1, 100]

/* ideal refresh rate */
int refreshRate = 60;

/* window content scale on x axis */
float xAxisWindowContentScale;
/* window content scale on y axis */
float yAxisWindowContentScale;

/* Main application window */
GLFWwindow *window;
/* Where data files live */
string RES_DIR = "";

/* shader for drawing with DIY colours */
shared_ptr<Program> progCol;
/* shader for picking */
shared_ptr<Program> progPick;
/* shader for heat solving */
shared_ptr<Program> progHeat;

// camera controls
/* arc ball */
ArcBall arcBall;
/* button state */
int buttonState = 0;
/* distance to object */
float distanceToObject = 12;
/* distance to object min */
const float minDistance = 1;
/* distance to object max */
const float maxDistance = 50;

/* VAO */
GLuint vao;
/* vbo position buffer ID */
GLuint posBufID;
/* location set in col_vert.glsl (or can be queried) */
GLuint aPosLocation = 0;
/* number of vertices */
const GLuint NumVertices = 3;
/* vertices buffer */
GLfloat vertices[NumVertices][3] = {
	{-1, -1, 0},
	{1, -1, 0},
	{1, 1, 0},
};

/* polygon soup */
shared_ptr<PolygonSoup> soup;
/* half edge data structure */
shared_ptr<HEDS> heds;
/* current half edge */
HalfEdge *currentHE;
/* which soup in list */
int whichSoup = 0;

/* number of soups */
const int numSoups = 16;
/* all soups, feel free to add your own, but control the size */
string soupFiles[numSoups] = {
	"bunnyLowRes.obj",
	"cow.obj",
	"cube.obj",
	"cube2obj.obj",
	"headtri.obj",
	"icosahedron.obj",
	"icoSphere2.obj",
	"icosphere6.obj",
	"sahedron.obj",
	"Sphere.obj",
	"torus.obj",
	"torusSmall.obj",
	"human.obj",			// quads
	"werewolf.obj",		// quads
	"bunny.obj",			// big challenge
	"better_dragon.obj",	// ultimate challenge
};

/* restart the heat */
bool restartHeatRequest = false;

/* click point in the window */
glm::vec2 clickPoint = glm::vec2(0, 0);
/* select heat source */
bool selectRequest = false;
/* selected heat source */
Vertex *selectedVertex = NULL;

/* mesh picking object */
MeshDrawPicking meshPicking;
/* mesh drawing object */
MeshDrawHeatGeo meshDraw;

/* new mesh loaded flag */
bool newMeshLoaded = false;

/* error callback */
static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

/* draw polygon soup as wire frame */
bool drawPolySoup = false;
/* draw half edge data structure */
bool drawHEDSMesh = true;
/* using wireframe mode */
bool wireFrame = false;
/* using face culling */
bool cullFace = true;
/* draw sample half edge */
bool drawHalfEdge = true;
/* run heat solver flag */
bool runHeatSolver = true;
/* step solver request */
bool stepHeatRequest = false;
/* draw picking flag */
bool drawPicking = false;
/* draw gradient flage */
bool drawGrad = false;

/* load a new mesh */
static void loadMesh(string filename)
{
	soup = NULL;
	soup = make_shared<PolygonSoup>(filename);
	heds = NULL;
	heds = make_shared<HEDS>(soup);

	heds->computeLaplacian();
	heds->precomputeQuantities();
	heds->initHeatFlow();

	newMeshLoaded = true;

	if (heds->faces->size() > 0) {
		currentHE = heds->faces->front()->he;
	} else {
		currentHE = NULL;
	}
}

/* key call back */
static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods) {
	if (action == GLFW_PRESS) {
		if (key == GLFW_KEY_ESCAPE) {
			glfwSetWindowShouldClose(window, GL_TRUE);
		} else if (key == GLFW_KEY_C) {
			cullFace = !cullFace;
		} else if (key == GLFW_KEY_SPACE) {
			if (currentHE->twin != NULL)
				currentHE = currentHE->twin;
		} else if (key == GLFW_KEY_N) {
			if (currentHE->next != NULL)
				currentHE = currentHE->next;
		} else if (key == GLFW_KEY_W) {
			if (currentHE->next != NULL)
				currentHE = currentHE->twin->next->twin->next->next;
		} else if (key == GLFW_KEY_PAGE_UP) {
			if (whichSoup > 0) {
				whichSoup--;
				loadMesh(RES_DIR + soupFiles[whichSoup]);
			} else {
				whichSoup = 0;
			}
		} else if (key == GLFW_KEY_PAGE_DOWN) {
			int soupListMax = numSoups - 1;
			if (whichSoup < soupListMax) {
				whichSoup++;
				loadMesh(RES_DIR + soupFiles[whichSoup]);
			} else {
				whichSoup = soupListMax;
			}
		} else if (key == GLFW_KEY_R) {
			restartHeatRequest = true;
		} else if (key == GLFW_KEY_S) {
			stepHeatRequest = true;
		} else if (key == GLFW_KEY_H) {
			drawHalfEdge = !drawHalfEdge;
		} else if (key == GLFW_KEY_F) {
			wireFrame = !wireFrame;
		} else if (key == GLFW_KEY_G) {
			drawGrad = !drawGrad;
		} else if (key == GLFW_KEY_P) {
			drawPicking = !drawPicking;
		} else if (key == GLFW_KEY_3) {
			drawPolySoup = !drawPolySoup;
		} else if (key == GLFW_KEY_M) {
			drawHEDSMesh = !drawHEDSMesh;
		} else if (key == GLFW_KEY_ENTER) {
			runHeatSolver = !runHeatSolver;
		} else if (key == GLFW_KEY_UP) {
			tSol *= 1.1;
			cout << "tsol = " << tSol << endl;
		} else if (key == GLFW_KEY_DOWN) {
			tSol /= 1.1;
			cout << "tsol = " << tSol << endl;
		}
	}
}

/* mouse event call back */
static void mouseDown_callback(GLFWwindow *window, int button, int action, int mods) {
	if (button == GLFW_MOUSE_BUTTON_LEFT) {
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) {
			buttonState = 0;
			return;
		}
		double x;
		double y;
		int windowWidth, windowHeight;
		glfwGetWindowSize(window, &windowWidth, &windowHeight);
		glfwGetCursorPos(window, &x, &y);
		arcBall.startRotation(x, y, windowWidth, windowHeight);
		buttonState = 1;
	} else if (button == GLFW_MOUSE_BUTTON_RIGHT) {
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_RELEASE) {
			double x;
			double y;
			glfwGetCursorPos(window, &x, &y);
			clickPoint.x = (float)x * xAxisWindowContentScale;
			clickPoint.y = (float)y * yAxisWindowContentScale;
			selectRequest = true;
		}
	}
}

/* mouse movement callback */
static void mouseMove_callback(GLFWwindow *window, double xpos, double ypos) {
	if (buttonState == 1) {
		double x;
		double y;
		int windowWidth, windowHeight;
		glfwGetWindowSize(window, &windowWidth, &windowHeight);
		glfwGetCursorPos(window, &x, &y);
		arcBall.updateRotation(x, y, windowWidth, windowHeight);
	}
}

/* mouse scrolling callback */
static void scroll_callback(GLFWwindow *window, double xoffset, double yoffset) {
	distanceToObject = glm::clamp((float)(distanceToObject - yoffset), minDistance, maxDistance);
}

/* put your initialization code here to improve performance! */
static void init() {
	GLSL::checkVersion();

	// Check how many texture units are supported in the vertex shader
	int tmp;
	glGetIntegerv(GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS, &tmp);
	cout << "GL_MAX_VERTEX_TEXTURE_IMAGE_UNITS = " << tmp << endl;
	// Check how many uniforms are supported in the vertex shader
	glGetIntegerv(GL_MAX_VERTEX_UNIFORM_COMPONENTS, &tmp);
	cout << "GL_MAX_VERTEX_UNIFORM_COMPONENTS = " << tmp << endl;
	glGetIntegerv(GL_MAX_VERTEX_ATTRIBS, &tmp);
	cout << "GL_MAX_VERTEX_ATTRIBS = " << tmp << endl;

	// Set background color.
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	// Enable z-buffer test.
	glEnable(GL_DEPTH_TEST);

	// setup OpenGL
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_LINE_SMOOTH);
	//glEnable(GL_POLYGON_SMOOTH);

	progCol = make_shared<Program>();
	progCol->setVerbose(true);
	progCol->setShaderNames(RES_DIR + "col_vert.glsl", RES_DIR + "col_frag.glsl");
	progCol->init();
	progCol->addUniform("P");
	progCol->addUniform("MV");
	progCol->addUniform("col");
	progCol->addAttribute("aPos");
	progCol->setVerbose(false);

	progPick = make_shared<Program>();
	progPick->setVerbose(true);
	progPick->setShaderNames(RES_DIR + "pick_vert.glsl", RES_DIR + "pick_frag.glsl");
	progPick->init();
	progPick->addUniform("P");
	progPick->addUniform("MV");
	progPick->addAttribute("vertex");
	progPick->addAttribute("color");
	progPick->setVerbose(false);

	progHeat = make_shared<Program>();
	progHeat->setVerbose(true);
	progHeat->setShaderNames(RES_DIR + "heat_vert.glsl", RES_DIR + "heat_frag.glsl");
	progHeat->init();
	progHeat->addUniform("P");
	progHeat->addUniform("MV");
	progHeat->addUniform("lightCamSpacePosition");
	progHeat->addUniform("lightColor");
	progHeat->addUniform("materialDiffuse");
	progHeat->addUniform("materialShininess");
	progHeat->addAttribute("vertex");
	progHeat->addAttribute("normal");
	progHeat->addAttribute("ut");
	progHeat->addAttribute("phi");
	progHeat->setVerbose(true);

	loadMesh(RES_DIR + soupFiles[0]);

	// mesh picking init
	meshPicking.init(progPick);

	// mesh draw init
	meshDraw.init(progHeat);

	// Create a buffers for doing some line drawing
	glGenVertexArrays(1, &vao);
	glBindVertexArray(vao);
	glGenBuffers(1, &posBufID);
	glBindBuffer(GL_ARRAY_BUFFER, posBufID);
	glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_DYNAMIC_DRAW);
	glEnableVertexAttribArray(aPosLocation);
	glVertexAttribPointer(aPosLocation, 3, GL_FLOAT, GL_FALSE, 0, (void *)0);

	// If there were any OpenGL errors, this will print something.
	// You can intersperse this line in your code to find the exact location
	// of your OpenGL error.
	GLSL::checkError(GET_FILE_LINE);
}

static void render() {
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	float aspect = width / (float)height;
	if (isnan(aspect)) {
		aspect = 0;
	}
	glViewport(0, 0, width, height);

	glLineWidth(1);

	// Clear framebuffer.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Create matrix stacks.
	shared_ptr<MatrixStack> P = make_shared<MatrixStack>();
	shared_ptr<MatrixStack> MV = make_shared<MatrixStack>();
	// Apply projection.
	P->pushMatrix();
	P->multMatrix(glm::perspective((float)(45.0 * M_PI / 180.0), aspect, 0.01f, 100.0f));
	// Apply camera transform.
	MV->pushMatrix();
	MV->translate(glm::vec3(0, 0, -distanceToObject));
	MV->multMatrix(arcBall.R);

	// if new mesh is loaded
	if (newMeshLoaded) {
		newMeshLoaded = false;
		meshPicking.setMesh(heds);
		meshDraw.setMesh(heds);
	}

	// use of face culling
	if (cullFace) {
		glEnable(GL_CULL_FACE);
	} else {
		glDisable(GL_CULL_FACE);
	}

	// draw in wireframe mode
	if (wireFrame) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); // set polygons to render in wire frame
	} else {
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL); // set polygons to render filled.
	}

	// restart heat simulation
	if (restartHeatRequest == true) {
		restartHeatRequest = false;
		heds->initHeatFlow();
	}

	// run heat solver
	if (runHeatSolver) {
		heds->solveHeatFlowStep(steps, tSol);
	} else if (stepHeatRequest) {
		stepHeatRequest = false;
		heds->solveHeatFlowStep(steps, tSol);
	}

	heds->updateGradu(); // gradient
	heds->updateDivx(); // divergence
	heds->solveDistanceStep(steps); // distance

	// select point
	glm::vec3 BBC;
	if (selectRequest)
	{
		selectRequest = false;
		int ID = meshPicking.pickTriangles(&clickPoint, &BBC, window, P, MV, progCol);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		cout << "ID = " << ID << " BBC = " << BBC.x << ", " << BBC.y << ", " << BBC.z << endl;
		if (ID >= 0) {
			// clear the constrained vertices
			for (auto &v : *(heds->vertices)) {
				v->constrained = false;
				v->u0 = 0;
			}
			// choose the closest vertex based on the barycentric coordinates
			if (BBC.x > BBC.y && BBC.x > BBC.z) {
				selectedVertex = heds->faces->at(ID)->he->head;
			} else if (BBC.y > BBC.x && BBC.y > BBC.z) {
				selectedVertex = heds->faces->at(ID)->he->next->head;
			} else {
				selectedVertex = heds->faces->at(ID)->he->next->next->head;
			}
			selectedVertex->constrained = true;
			selectedVertex->u0 = 1;
			selectedVertex->ut = 1;
		}
	}

	// draw HEDS mesh
	if (drawHEDSMesh) {
		if (drawPicking) {
			meshPicking.pickTriangles(&clickPoint, &BBC, window, P, MV, progCol);
		} else {
			meshDraw.drawVBOs(P, MV);
		}
	}

	// draw soup wireframe
	if (drawPolySoup) {
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		soup->display(progCol, P, MV);
	}

	// draw sample halfedge
	if (drawHalfEdge) {
		if (currentHE != NULL) {
			currentHE->display(progCol, P, MV, vao, posBufID);
		}
	}

	// draw gradients
	if (drawGrad) {
		meshDraw.drawGrad(P, MV, progCol);
	}

	// Pop matrix stacks.
	MV->popMatrix();
	P->popMatrix();

	GLSL::checkError(GET_FILE_LINE);
}

int main(int argc, char **argv) {
	if (argc < 2) {
		cout << "Please specify the resource directory." << endl;
		return 0;
	}
	RES_DIR = argv[1] + string("/");

	// Set error callback.
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if (!glfwInit()) {
		return -1;
	}

	// https://en.wikipedia.org/wiki/OpenGL
	// hint to use OpenGL 4.1 on all paltforms
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(640, 480, "LEOPOLDO ZUGASTI 260919951", NULL, NULL);
	if (!window) {
		glfwTerminate();
		return -1;
	}
    // get the window content scales
    glfwGetWindowContentScale(window, &xAxisWindowContentScale, &yAxisWindowContentScale);
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if (glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}
	glGetError(); // A bug in glewInit() causes an error that we can safely ignore.
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
	// Set vsync.
	GLFWmonitor *monitor = glfwGetPrimaryMonitor();
	const GLFWvidmode *mode = glfwGetVideoMode(monitor);
	refreshRate = mode->refreshRate;
	glfwSwapInterval((int)((float)(refreshRate) / 60.0));
	// Set keyboard callback.
	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mouseDown_callback);
	glfwSetCursorPosCallback(window, mouseMove_callback);
	glfwSetScrollCallback(window, scroll_callback);
	// Initialize scene.
	init();
	// Loop until the user closes the window.
	while (!glfwWindowShouldClose(window)) {
		// Render scene.
		render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();
	}
	// Quit program.
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
