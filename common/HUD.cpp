////////////////////////////////////////////////////////////////////////////////
// A quick-and-dirty Heads-Up-Display class for GLUT applications.
// 
// This is a pretty minimal HUD class for openGL applications that use the GLUT
// window manager. The HUD contains an array of HUDElements as a private data member,
// and a private integer that keeps track of the last added element. New elements
// are added with HUD::AddElement(). The string, position, and color of existing
// elements are set with the public HUD::set*() methods. A call to HUD::Render()
// prints all the elements (most of which are normally empty strings) in the current
// GLUT window. A call to HUD::Clear() will set the string value of all HUDElements
// to "".
// 
// Usage:
// HUD hud;
// int e1 = hud.AddElement("Mitzi is Pooch",0.5,0.5); // add new element
// hud.SetElement(e1,"She is so Pooch!"); // modify existing element
// 
// Notes:
// This is a bare bones implementations. There are almost no checks, and certainly
// no error messages. Some potential gotchas to keep in mind are:
//   1) The string manipulation is done with char pointers and the unsecure strcpy().
//      When using the Microsoft VC++ compiler you may #define _CRT_SECURE_NO_WARNINGS
//      to suppress security warnings. It is the user's responsibility to never send
//      the HUD::Set* functions anything but a valid pointer to a null terminated string,
//      no longer than HUD_ELEMENT_LENGTH.
//   2) All HUDElements strings in HUD::elements[] are implicitly initialized to
//      an empty string. The HUD::Render() function renders ALL elements, not just
//      ones you previously set.
//   3) There is no HUD::Delete() method. You may set an element's string to ""
//      and keep adding new elements as needed.
//   4) Calling HUD::AddElement() more than HUD_MAX_ELEMENTS times in a program
//      has no effect.
//   5) Calling HUD::Set*() with an index greater than HUD::LAST_ELEMENT has
//      no effect.
// 
// Author: Naor Movshovitz
// Edits: $Date: 2012-09-19 12:56:24 -0700 (Wed, 19 Sep 2012) $ $Author: nmovshov $
///////////////////////////////////////////////////////////////////////////////

#include "HUD.h"

// Request lint level warnings with the LINT macro on Microsoft compilers
#ifdef LINT
#pragma warning(push,4)
#pragma warning(disable:4100) // unreferenced formal parameter
#endif

HUD::HUD()
{
	LAST_ELEMENT = 0;
}
unsigned int HUD::AddElement(const char* s, const float xpos, const float ypos)
{
	LAST_ELEMENT++;
	if (LAST_ELEMENT<HUD_MAX_ELEMENTS)
		SetElement(LAST_ELEMENT,s,xpos,ypos);
	return LAST_ELEMENT;
}
void HUD::SetElement(const unsigned int ind, const char* s, const float xpos, const float ypos)
{
	if (ind<=LAST_ELEMENT)
	{
		strcpy(elements[ind].m_string,s);
		elements[ind].xpos = xpos * glutGet(GLUT_WINDOW_WIDTH);
		elements[ind].ypos = ypos * glutGet(GLUT_WINDOW_HEIGHT);
	}
}
void HUD::SetElement(const unsigned int ind, const char* s)
{
	if (ind<=LAST_ELEMENT)
	{
		strcpy(elements[ind].m_string,s);
	}
}
void HUD::SetElementColor(const unsigned int ind, const float r, const float g, const float b, const float alpha/*=1.0*/)
{
	if (ind<=LAST_ELEMENT)
	{
		elements[ind].red   = r;
		elements[ind].green = g;
		elements[ind].blue  = b;
		elements[ind].alpha = alpha;
	}
}
void HUD::Render()
{
	// First, switch to orthographic projections
	glMatrixMode(GL_PROJECTION);	// Switch to projection mode
	glPushMatrix();					// Save previous matrix which contains the settings for the perspective projection
	glLoadIdentity();				// Reset matrix
	gluOrtho2D(0, glutGet(GLUT_WINDOW_WIDTH), glutGet(GLUT_WINDOW_HEIGHT), 0);	// Set a 2D orthographic projection
	glDisable(GL_LIGHTING);			// HUD elements are not subject to lighting
	glMatrixMode(GL_MODELVIEW);		// Switch back to model view mode

	// Then render all HUD elements with bitmap fonts
	for (unsigned int k=0; k<HUD_MAX_ELEMENTS; k++)
	{
		glPushMatrix();
		glColor4f(elements[k].red,elements[k].green,elements[k].blue,elements[k].alpha);
		glLoadIdentity();
		glRasterPos2f(elements[k].xpos,elements[k].ypos);
		for (char *c=elements[k].m_string; *c!='\0'; c++)
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18,*c);
		glPopMatrix();
	}
	
	// And restore perspective projection
	glMatrixMode(GL_PROJECTION);
	glPopMatrix(); // Restore previously saved projection settings
	glMatrixMode(GL_MODELVIEW);
	glEnable(GL_LIGHTING); // (don't forget to turn ON the lights when you leave)
}
void HUD::Clear()
{
	int le=LAST_ELEMENT;
	LAST_ELEMENT = HUD_MAX_ELEMENTS;
	for (unsigned int k=0; k<HUD_MAX_ELEMENTS; k++) SetElement(k,"",0.0f,0.0f);
	LAST_ELEMENT=le;
}

// End lint level warnings
#ifdef LINT
#pragma warning(pop)
#endif