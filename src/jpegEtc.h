#ifdef _WIN32

#pragma once

// deliberately does not use util.h
// because needs all Windows in the implementation
//
#include <string>
#include "util_nocopy.h"

class JPEGEtcGrey : public NoCopy {
public:
	JPEGEtcGrey() {}
	virtual ~JPEGEtcGrey() {}
	virtual float intensity(int w, int h) const = 0;
};

class JPEGEtcRGB : public NoCopy {
public:
	JPEGEtcRGB() {}
	struct Color {
		unsigned char r, g, b;  
		Color() {}
		Color(unsigned char r, unsigned char g, unsigned char b) : r(r), g(g), b(b) {}
	};
	virtual ~JPEGEtcRGB() {}
	virtual Color rgb(int w, int h) const = 0;
};

class JPEGEtc : public NoCopy {
public:
	static JPEGEtc* alloc  ();
	static void     dealloc(JPEGEtc*& made);

	JPEGEtc() {}
	virtual ~JPEGEtc() {}
	virtual void listEncoders() = 0;

	class Photo : public NoCopy {
	public: 
		Photo() {}
		virtual ~Photo() {}
		virtual void saveAsBMP(std::string fileName) = 0;

		virtual void setAllPixels     (JPEGEtcGrey const & input) = 0;
		virtual void setNonzeroPixels (JPEGEtcRGB const & input) = 0;
		virtual void setOneGreyPixel  (float intensity01, int w, int h) = 0;
	};

	virtual Photo* makePhoto(size_t width, size_t height) = 0;
	virtual void   freePhoto(Photo*& made) = 0;
};

void jpegEtcTest();

#endif