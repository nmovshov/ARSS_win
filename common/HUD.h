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
// This is a bare bones implementation. There are almost no checks, and certainly
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
// Edits: $Date: 2013-04-06 20:11:42 -0700 (Sat, 06 Apr 2013) $ $Author: nmovshov $
///////////////////////////////////////////////////////////////////////////////

#ifndef HUD_H
#define HUD_H

#ifdef __APPLE__
#include "GLUT/glut.h"
#else
#include "GL/glut.h"					// The GLUT window and event manager
#endif
#include "ncclib.h"

#define HUD_MAX_ELEMENTS 64
#define HUD_ELEMENT_LENGTH 256

class HUD;

class HUDElement
{
private:
	char*	m_string;
	float	xpos, ypos;
	float	red, blue, green, alpha;

	HUDElement() {
		m_string = new char[HUD_ELEMENT_LENGTH];
		strcpy(m_string,"");
		xpos=ypos=0.0f;
		red=green=blue=alpha=1.0f;
	}
public:
	~HUDElement() {
		delete [] m_string;
	}

	friend class HUD;
};

class HUD
{
private:
	HUDElement elements[HUD_MAX_ELEMENTS];
	unsigned int LAST_ELEMENT;
public:
	HUD();
	unsigned int AddElement(const char* s, const float xpos, const float ypos);
	void SetElement(const unsigned int ind, const char* s, const float xpos, const float ypos);
	void SetElement(const unsigned int ind, const char* s);
	void SetElementColor(const unsigned int ind, const float r, const float g, const float b, const float alpha=1.0);
	void Clear();
	void Render();
};

#endif