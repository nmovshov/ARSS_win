///////////////////////////////////////////////////////////////////////////////
// This file contains the hi-, mid-, and low-level routines used to draw things
// in the ARSS solution. Only the high-level routines are exposed (in the header)
// and allowed to be externally linked. The high-level routines are:
//	- DrawShapes(PxShape*) - accepts a pointer to a PxShape object, sets some 
//	  attributes (i.e. color) based on the object's userData, and dispatches the
//	  object to specilized Draw*() functions based on PxShape::getGeometryType().
//	- Draw*Grid - routines to draw various grids for visual debugging.
//	- DrawLine.
//	- DrawArrow.
//	
//	The mid-level functions are:
//	- Draw*Shape(PxShape*). These functions use shape data to construct the necessary
//	  OpenGL transformation and then call a low-level routine to issue the OpenGL
//	  drawing commands.
//	  
//	The low-level functions are:
//	- Render*().
//
// Author: Naor
// Changed: $Date: 2013-02-19 19:50:50 -0800 (Tue, 19 Feb 2013) $ $Author: nmovshov $
///////////////////////////////////////////////////////////////////////////////

#include "Drawers.h"

/* I use forward declarations for all mid- and low-level routines, instead of a header,
 * to avoid the temptation of calling these functions directly from the application.
 * The calls wouldn't work anyway because I suppress external linkage. This allows
 * me some leeway to use "unsafe" techniques in this file.
 */
static void DrawBox(PxShape*);
static void DrawSphere(PxShape*);
static void DrawPlane(PxShape*);
static void DrawCapsule(PxShape*);
static void DrawConvex(PxShape*);
static void RenderBox();
static void RenderSphere();
static void RenderCylinder();
static void RenderPlane();
static void RenderRightHandedAxes();

// Top level "dispatcher" function available for external linkage
void DrawShape(PxShape* shape)
{
	// Start by clearing the top of the matrix stack so that mid level routine can count on a clean slate
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();

	// Also save the current attributes so we can mess with colors etc.
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);

	// A vector<GLUbyte> in the shape's userData overrides the current (actor) color
	if (shape->userData) {
		glColor3ubv(static_cast<GLubyte*>(shape->userData));
	}

	// Next, dispatch the shape to an appropriate drawing routine
	PxGeometryHolder geom = shape->getGeometry();
	switch (geom.getType())
	{
	case PxGeometryType::eBOX:
		DrawBox(shape);
		break;
	case PxGeometryType::eSPHERE:
		DrawSphere(shape);
		break;
	case PxGeometryType::eCAPSULE:
		DrawCapsule(shape);
		break;
	case PxGeometryType::ePLANE:
		DrawPlane(shape);
		break;
	case PxGeometryType::eCONVEXMESH:
		DrawConvex(shape);
		break;
	default:
		printf("Unknown shape cannot be rendered!\a");
		exit(1);
		break;
	}

	// Finally, return the stacks to whatever they contained in the main application (probably nothing)
	glPopMatrix();
	glPopAttrib();
}

// Other hi-level, exposed, externally linked routines
void Draw3DBoxGrid()
{
	// Grid lines coordinates
	static const GLfloat verts[] = {-1.0f,+0.0f,+0.0f,
									+1.0f,+0.0f,+0.0f,
									+0.0f,-1.0f,+0.0f,
									+0.0f,+1.0f,+0.0f,
									+0.0f,+0.0f,-1.0f,
									+0.0f,+0.0f,+1.0f,
									};
	glVertexPointer(3,GL_FLOAT,0,&(verts[0]));

	// Save and clear the stacks
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);

	// Draw the grid
	glDisable(GL_LIGHTING);
	glColor3ubv(ncc::rgb::gChartreuse);
	for (int j=-10; j<11; j++) {
		for (int k=-10; k<11; k++) {
			for (int l=-10; l<11; l++) {
				glTranslatef(j,k+0.01,l);
				glDrawArrays(GL_LINES,0,6);
				glLoadIdentity();
			}
		}
	}

	// Restore the stacks
	glPopMatrix();
	glPopAttrib();
}
void DrawXZGrid()
{
	// Save and clear the stacks
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);

	// Draw the grid
	glDisable(GL_LIGHTING);
	glColor3ubv(ncc::rgb::gChartreuse);
	glTranslated(0.0,0.01,0.0); // line "painted on" the ground
	glBegin(GL_LINES);
	for (int k=-10; k<11; k++)
	{
		glVertex3i(k,0,-10);
		glVertex3i(k,0,+10);
		glVertex3i(-10,0,k);
		glVertex3i(+10,0,k);
	}
	glEnd();

	// Restore the stacks
	glPopMatrix();
	glPopAttrib();
}
void DrawXYGrid()
{
	// Save and clear the stacks
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT);

	// Draw the grid
	glDisable(GL_LIGHTING);
	glColor3ubv(ncc::rgb::gChartreuse);
	glTranslated(0.0,0.0,0.01); // line "painted on" the XY plane
	glBegin(GL_LINES);
	for (int k=-10; k<11; k++)
	{
		glVertex3i(k,-10,0);
		glVertex3i(k,+10,0);
		glVertex3i(-10,k,0);
		glVertex3i(+10,k,0);
	}
	glEnd();

	// Restore the stacks
	glPopMatrix();
	glPopAttrib();
}
void DrawActorAxes( PxRigidActor* actor, GLfloat scale/*=1.0f*/ )
{
	// Save and clear the stacks
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_LINE_BIT);

	// Prepare
	PxTransform pose;
	PxReal angle;
	PxVec3 axis;

	// Inquire
	pose = actor->getGlobalPose();
	pose.q.toRadiansAndUnitAxis(angle,axis);
	angle=angle*180.0f/PxPi;

	// Translate, rotate, scale (in that order!!!!)
	glTranslatef(pose.p.x,pose.p.y,pose.p.z);
	glRotatef(angle,axis.x,axis.y,axis.z);
	glScalef(scale,scale,scale);

	// Render
	glDisable(GL_LIGHTING);
	RenderRightHandedAxes();

	// Restore the stacks
	glPopMatrix();
	glPopAttrib();
}
void DrawLine(PxVec3 p1, PxVec3 p2, GLfloat thickness/*=1*/, const GLubyte* color/*=NULL*/)
{
	// Save and clear the stacks
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_LINE_BIT);

	// Set up properties
	glDisable(GL_LIGHTING);
	glLineWidth(thickness);
	if (color)
		glColor3ubv(color);

	// Draw the line
	glBegin(GL_LINES);
		glVertex3f(p1.x,p1.y,p1.z);
		glVertex3f(p2.x,p2.y,p2.z);
	glEnd();

	// Restore the stacks
	glPopMatrix();
	glPopAttrib();
}
void DrawArrow(PxVec3 base, PxVec3 head, GLfloat headSize/*=0.12*/, GLfloat thickness/*=4*/, const GLubyte* color/*=NULL*/)
{
	// Save and clear the stacks
	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glLoadIdentity();
	glPushAttrib(GL_CURRENT_BIT | GL_ENABLE_BIT | GL_LIGHTING_BIT | GL_LINE_BIT);

	// Set up properties
	glDisable(GL_LIGHTING);
	glLineWidth(thickness);
	if (color)
		glColor3ubv(color);

	// Draw the shaft
	glBegin(GL_LINES);
		glVertex3f(base.x,base.y,base.z);
		glVertex3f(head.x,head.y,head.z);
	glEnd();
	
	// Draw the head
	PxVec3 n = (head-base); n.normalize();
	PxVec3 z(0.0f,0.0f,1.0f);
	PxVec3 axis = z.cross(n);
	PxReal angle = PxAcos(n.dot(z)) * 180.0f / PxPi;
	headSize = headSize * (head-base).magnitude();
	glTranslatef(head.x,head.y,head.z);
	if (angle)
		glRotatef(angle,axis.x,axis.y,axis.z);
	glutSolidCone(0.5*headSize,headSize,12,12);

	// Restore the stacks
	glPopMatrix();
	glPopAttrib();
}


// Mid level model and view transformation functions
static void DrawBox(PxShape* shape)
{
	// Prepare
	PxBoxGeometry	geom;
	PxTransform		pose;
	PxReal			angle;
	PxVec3			axis;

	// Inquire
	shape->getBoxGeometry(geom);
	pose = PxShapeExt::getGlobalPose(*shape);
	pose.q.toRadiansAndUnitAxis(angle,axis);
	angle=angle*180.0f/PxPi;

	// Translate, rotate, scale (in that order!!!!)
	glTranslatef(pose.p.x,pose.p.y,pose.p.z);
	glRotatef(angle,axis.x,axis.y,axis.z);
	glScalef(geom.halfExtents.x,geom.halfExtents.y,geom.halfExtents.z);

	// Draw
	RenderBox();
}
static void DrawSphere(PxShape* shape)
{
	// Prepare
	PxSphereGeometry	geom;
	PxTransform			pose;
	
	// Inquire
	shape->getSphereGeometry(geom);
	pose = PxShapeExt::getGlobalPose(*shape);

	// Translate, scale (in that order!!!!)
	glTranslatef(pose.p.x,pose.p.y,pose.p.z);
	glScalef(geom.radius,geom.radius,geom.radius);

	// Draw
	RenderSphere();
}
static void DrawCapsule(PxShape* shape)
{
	// Prepare
	PxCapsuleGeometry	geom;
	PxTransform			pose;
	PxReal				angle;
	PxVec3				axis;
	GLUquadricObj*		qobj;

	// Inquire
	shape->getCapsuleGeometry(geom);
	pose = PxShapeExt::getGlobalPose(*shape);
	pose.q.toRadiansAndUnitAxis(angle,axis);
	angle=angle*180.0f/PxPi;

	// Translate, rotate (in that order!!!!)
	glTranslatef(pose.p.x,pose.p.y,pose.p.z);
	glRotatef(angle,axis.x,axis.y,axis.z);

	// Draw top cap
	glPushMatrix();
	glTranslatef(+geom.halfHeight,0.0f,0.0f);
	glScalef(geom.radius,geom.radius,geom.radius);
	RenderSphere();
	glPopMatrix();
	// Draw bottom cap
	glPushMatrix();
	glTranslatef(-geom.halfHeight,0.0f,0.0f);
	glScalef(geom.radius,geom.radius,geom.radius);
	RenderSphere();
	glPopMatrix();
	// Draw cylinder
	qobj = gluNewQuadric();
	gluQuadricDrawStyle(qobj,GLU_FILL);
	gluQuadricNormals(qobj,GLU_SMOOTH);
	glPushMatrix();
	glTranslatef(-geom.halfHeight,0,0);
	glRotatef(90,0,1,0);
	gluCylinder(qobj,geom.radius,geom.radius,2.02f*geom.halfHeight,16,6);
	glPopMatrix();
	gluDeleteQuadric(qobj);
}
static void DrawConvex(PxShape* shape)
{
	// Prepare
	PxConvexMeshGeometry	geom;
	PxConvexMesh*			mesh;
	PxTransform				pose;
	PxReal					angle;
	PxVec3					axis;

	// Inquire
	shape->getConvexMeshGeometry(geom);
	mesh = geom.convexMesh;
	pose = PxShapeExt::getGlobalPose(*shape);
	pose.q.toRadiansAndUnitAxis(angle,axis);
	angle=angle*180.0f/PxPi;

	// Translate, rotate, scale (in that order!!!!)
	glTranslatef(pose.p.x,pose.p.y,pose.p.z);
	glRotatef(angle,axis.x,axis.y,axis.z);
	glScalef(geom.scale.scale.x,geom.scale.scale.y,geom.scale.scale.z);

	// Iterate over faces to triangulate and draw
	PxU32			nbConvexVerts = mesh->getNbVertices();
	PxU32			nbConvexPolys = mesh->getNbPolygons();
	const PxVec3*	convexVerts	= mesh->getVertices();
	const PxU8*		indexBuffer = mesh->getIndexBuffer(); // This is a long list of indices for face vertices
	
	for (PxU32 k=0; k<nbConvexPolys; k++)
	{
		PxHullPolygon face;
		if (mesh->getPolygonData(k,face))
		{
			glBegin(GL_POLYGON);
			glNormal3f(face.mPlane[0],face.mPlane[1],face.mPlane[2]);
			const PxU8* faceIdx = indexBuffer + face.mIndexBase;
			for (PxU32 j=0; j<face.mNbVerts; j++)
			{
				PxVec3 v = convexVerts[faceIdx[j]];
				glVertex3f(v.x,v.y,v.z);
			}
			glEnd();
		}
	}
}
static void DrawPlane(PxShape* shape)
{
	// Prepare
	PxPlaneGeometry	geom;
	PxTransform		pose;
	PxReal			angle;
	PxVec3			axis;

	// Inquire
	shape->getPlaneGeometry(geom);
	pose = PxShapeExt::getGlobalPose(*shape);
	pose.q.toRadiansAndUnitAxis(angle,axis);
	angle=angle*180.0f/PxPi;

	// Translate, rotate, scale (in that order!!!)
	glTranslatef(pose.p.x,pose.p.y,pose.p.z);
	glRotatef(angle,axis.x,axis.y,axis.z);
	glScalef(1024.0f,1024.0f,1024.0f);

	// Draw
	RenderPlane();
}

// Low level vertex assignment functions
static void RenderBox()
{
	glutSolidCube(2.0);
}
static void RenderSphere()
{
	glutSolidSphere(1.0,12,12);
}
static void RenderPlane()
{
	// a "neutral" PxPlane geometry is the YZ plane pointing right
	static const GLfloat verts[] =	{	0.0f,-1.0f,-1.0f,
										0.0f,-1.0f,+1.0f,
										0.0f,+1.0f,+1.0f,
										0.0f,-1.0f,-1.0f,
										0.0f,+1.0f,+1.0f,
										0.0f,+1.0f,-1.0f
									};
	glNormal3f(1.0f,0.0f,0.0f);
	glVertexPointer(3,GL_FLOAT,0,verts);
	glDisable(GL_LIGHTING);
	glDrawArrays(GL_TRIANGLES,0,6);
	glEnable(GL_LIGHTING);
}
static void RenderRightHandedAxes()
{
	glBegin(GL_LINES);
		// X-axis
		glColor3f(1.0f,0.0f,0.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(1.0f,0.0f,0.0f);

		// Y-Axis
		glColor3f(0.0f,1.0f,0.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.01,1.0f,0.0f);

		// Z-axis
		glColor3f(0.0f,0.0f,1.0f);
		glVertex3f(0.0f,0.0f,0.0f);
		glVertex3f(0.0f,0.0f,1.0f);
	glEnd();
}