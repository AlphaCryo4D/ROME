#ifdef _WIN32

// GDI - https://msdn.microsoft.com/en-us/library/vs/alm/ms533842(v=vs.85).aspx
//
// Need all Windows support, not just WIN32_LEAN_AND_MEAN as defined in util.h
//
#ifdef _WIN32
#define JPEGETC_IMPLEMENTED 
#endif

#ifdef JPEGETC_IMPLEMENTED
#undef NOMINMAX
#undef WIN32_LEAN_AND_MEAN
#include <Windows.h>
#include <gdiplus.h>
#undef max
#undef min
#endif

#include "./jpegEtc.h"

#include <algorithm>
#include <iostream>
#include <string>
#include <assert.h>

inline float square(float x) { 
	auto result = x*x;
	return result; 
}

#if defined(JPEGETC_IMPLEMENTED) && defined(_WIN32)

using namespace Gdiplus;
#pragma comment (lib,"Gdiplus.lib")

// https ://msdn.microsoft.com/en-us/library/ms235631.aspx
class WString {
public:
	size_t   const buflen;
	wchar_t* const bufptr;

	WString(std::string const & rhs) : buflen(rhs.size()+1), bufptr(new wchar_t[buflen]) {
		size_t convertedChars = 0;
		mbstowcs_s(&convertedChars, bufptr, buflen, rhs.c_str(), _TRUNCATE);
	}

	~WString() { delete bufptr; }
};
#endif


#if defined(JPEGETC_IMPLEMENTED) && defined(_WIN32)
class JPEGEtcImpl : public JPEGEtc {
	GdiplusStartupInput gdiplusStartupInput;
	ULONG_PTR			gdiplusToken;
public:
	JPEGEtcImpl() : JPEGEtc(), jpgClsidFound(false) {
		// Initialize GDI+.
		GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, NULL);
	}
	virtual ~JPEGEtcImpl() {
		GdiplusShutdown(gdiplusToken);
	}

	class ScanEncodersCallable {
	public:
		virtual ~ScanEncodersCallable() {}
		virtual bool operator()(ImageCodecInfo& ici) {
			return false;
		}
	};
	void scanEncoders(ScanEncodersCallable& callable) {
		UINT  num;
		UINT  size;
		GetImageEncodersSize(&num, &size);
		char* p = (new char[size]);
		auto pImageCodecInfo = (ImageCodecInfo*)p;
		GetImageEncoders(num, size, pImageCodecInfo);
		for (UINT j = 0; j < num; ++j) {
			ImageCodecInfo& ici = pImageCodecInfo[j];
			if (!callable(ici)) break;
		}
		pImageCodecInfo = nullptr;
		delete [] p;
	}

	void listEncoders() {
		class Callable : public ScanEncodersCallable {
		public:
			virtual bool operator()(ImageCodecInfo& ici) {
				wprintf(L"%s\n", ici.MimeType);
				return true;
			}
		} callable;
		scanEncoders(callable);
	}

	CLSID jpgClsid;
	bool  jpgClsidFound;

	class PhotoImpl : public Photo {
		JPEGEtcImpl* const parent;
		Gdiplus::Bitmap gdiImage;

	public:
		static PhotoImpl* alloc(JPEGEtcImpl* const parent, size_t width, size_t height) {
			return new PhotoImpl(parent, width, height);
		}
		PhotoImpl(JPEGEtcImpl* const parent, size_t width, size_t height)
			: Photo(), parent(parent), gdiImage(INT(width), INT(height), PixelFormat16bppRGB555)
		{}
		void saveAsBMP(std::string fileName) {

			fileName += ".bmp";

			// https://msdn.microsoft.com/en-us/library/windows/desktop/ms535407(v=vs.85).aspx
			WString wfileName(fileName);

			if (!parent->jpgClsidFound) {
				class Callable : public ScanEncodersCallable {
					JPEGEtcImpl* const parent;
				public:
					Callable(JPEGEtcImpl* const parent) : parent(parent) {}
					CLSID pngClsid;
					virtual bool operator()(ImageCodecInfo& ici) {
						auto result = wcscmp(ici.MimeType, L"image/bmp");
						//std::wcerr << L"Consider " << ici.MimeType << " result:" << result << std::endl;
						if (result != 0) {
							//std::wcerr << L"    no match" << std::endl;
							return true;
						}
						//std::wcerr << L"    match!!!" << std::endl;
						parent->jpgClsid = ici.Clsid;
						parent->jpgClsidFound = true;
						return false;
					}
				} callable(parent);
				parent->scanEncoders(callable);
				if (!parent->jpgClsidFound) {
					std::cerr << "Could not find image/bmp encoder to encode" << fileName << std::endl;
					return;
				}
				else {
					std::cout << "Found codec" << std::endl;
				}
			}

			Gdiplus::Status status = gdiImage.Save(wfileName.bufptr, &parent->jpgClsid);
			std::cerr << "Save as " << fileName << " status:" << status << std::endl;
		}

		virtual void setNonzeroPixels(JPEGEtcRGB const & input) {
			UINT const width  = gdiImage.GetWidth();
			UINT const height = gdiImage.GetHeight();
			for (UINT w = 0; w < width; w++) {
				for (UINT h = 0; h < height; h++) {
					auto rgb = input.rgb(int(w),int(h));
					auto b   = rgb.b;
					auto g   = rgb.g;
					auto r   = rgb.r;
					if (r + b + g == 0) continue;
					Gdiplus::Color color(r,g,b);
					Gdiplus::Status status = gdiImage.SetPixel(w, h, color);
					if (status != 0) {
						std::cerr << "gdiImage.SetPixel returns " << status << std::endl;
						return;
					}
				}
			}
		}

		virtual void setOneGreyPixel(float intensity01, int w, int h) {
			assert(0.0f <= intensity01);
			assert(1.0f >= intensity01);
			BYTE d = BYTE(255) - BYTE(intensity01 * 255);
			// std::cerr << int(d) << ' ';
			Gdiplus::Color color(d, d, d);
			Gdiplus::Status status = gdiImage.SetPixel(w, h, color);
			if (status != 0) {
				std::cerr << "gdiImage.SetPixel returns " << status
					<< " w:" << w << " h:" << h << " d:" << int(d) << std::endl;
				return;
			}
		}

		virtual void setAllPixels(JPEGEtcGrey const & input) {
			INT const width  = gdiImage.GetWidth();
			INT const height = gdiImage.GetHeight();
			// std::cerr << "setAllPixels width:" << width << " height:" << height << std::endl;
			auto fMin = input.intensity(0, 0); auto fMax = fMin;
			for (INT w = 0; w < width; w++) {
				for (INT h = 0; h < height; h++) {
					auto f = input.intensity(int(w), int(h));
					if (f < 0) f = 0;
					fMin = std::min(fMin, f);
					fMax = std::max(fMax, f);
				}
			}
			// std::cerr << "Encoded fMin:" << fMin << " fMax:" << fMax << std::endl;
			if (fMax == fMin) fMax = fMin + 1;
			for (INT w = 0; w < width; w++) {
				for (INT h = 0; h < height; h++) {
					auto f = input.intensity(int(w), int(h));
					f = (f - fMin) / (fMax - fMin);
					setOneGreyPixel(f, w, h);
				}
				// std::cerr << std::endl;
			}
			// std::cerr << "Encoded " << width*height << " pixels" << std::endl;
		}
	};

	virtual Photo* makePhoto(size_t width, size_t height) { return PhotoImpl::alloc(this, width, height); }
	virtual void   freePhoto(Photo*& made) { delete made; made = nullptr;  }
};
#endif

JPEGEtc* JPEGEtc::alloc() { 
	return 
#if defined(JPEGETC_IMPLEMENTED) && defined(_WIN32)
		new JPEGEtcImpl;
#else
		nullptr;
#endif
}
void     JPEGEtc::dealloc(JPEGEtc*& made) { delete made; made = nullptr; }


class JPEGEtcTest_DrawPhoto : public JPEGEtcGrey {
public:
	JPEGEtcTest_DrawPhoto(size_t width, size_t height) : JPEGEtcGrey(), width(width), height(height) {
		const int h = int(height);
		const int w = int(width);
		maxIntensity = 1.0;
		maxIntensity = 
			std::max(std::max(intensity(0,  0),intensity(w-1,  0)),
				     std::max(intensity(0,h-1),intensity(w-1,h-1)));
	}
	virtual float intensity(int w, int h) const {
		auto result = (square(w-int(width /2))+
			           square(h-int(height/2))
			          )/maxIntensity;
		return result;
	}
private:
	const size_t width, height;
	float  maxIntensity;
};

void jpegEtcTest() {
	const char* bmpFilename = "C:/TEMP/rome_map32_jpegEtcTest";
	const size_t w = 100, h = 100;
	if (auto jpegEtc = JPEGEtc::alloc()) {
		std::cerr << __FILE__ << " jpegEtcTest writing " << bmpFilename << std::endl;
		auto photo = jpegEtc->makePhoto(w,h);
		JPEGEtcTest_DrawPhoto brush(w, h);
		photo->setAllPixels(brush);
		photo->saveAsBMP(bmpFilename);
		jpegEtc->freePhoto(photo);
		JPEGEtc::dealloc(jpegEtc);
	} else {
		std::cerr << __FILE__ << " jpegEtc NYI" << std::endl;
	}
}

#endif