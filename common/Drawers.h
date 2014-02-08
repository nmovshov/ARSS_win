///////////////////////////////////////////////////////////////////////////////
// This file headers the hi-level routines used to draw things in the ARSS solution.
// Only the high-level routines are exposed here and allowed to be externally linked.
// The high-level routines are:
//	- DrawShapes(PxShape*) - accepts a pointer to a PxShape object, sets some 
//	  attributes (i.e. color) based on the object's userData, and dispatches the
//	  object to specilized Draw*() functions based on PxShape::getGeometryType().
//	- Draw*Grid - Routines to draw various grids for visual debugging.
//	- DrawLine.
//	- DrawArrow.
//
// Author: Naor
// Changed: $Date: 2012-10-01 11:59:37 -0700 (Mon, 01 Oct 2012) $ $Author: nmovshov $
///////////////////////////////////////////////////////////////////////////////

#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include "GL/glut.h"					// The GLUT window and event manager
#endif
#include "ncclib.h"
#include "PxPhysicsAPI.h"
using namespace physx;

void DrawShape(PxShape*);
void Draw3DBoxGrid();
void DrawXZGrid();
void DrawXYGrid();
void DrawLine(PxVec3 p1, PxVec3 p2, GLfloat thickness=1, const GLubyte* color=NULL);
void DrawArrow(PxVec3 base, PxVec3 head, GLfloat headSize=0.12, GLfloat thickness=4, const GLubyte* color=NULL);
void DrawActorAxes(PxRigidActor* actor, GLfloat scale=1.0f);